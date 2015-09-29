#!/usr/bin/env python
from __future__ import print_function

import os
from platform import system
import sys

import common
from tasks import get_blast_format_task, \
                  get_blast_task, \
                  get_transcriptome_stats_task, \
                  get_busco_task, \
                  get_group_task, \
                  get_link_file_task, \
                  get_transdecoder_predict_task, \
                  get_transdecoder_orf_task, \
                  get_hmmscan_task, \
                  get_cmscan_task, \
                  get_lastal_task

def get_annotate_tasks(transcriptome, prog_paths, database_dict, 
                       n_threads=1, user_databases=[]):

    tasks = []
    results = {}

    user_database_dict = {}
    if user_databases:
        for db in user_databases:
            db_path = db + '.db'
            tasks.append(
                get_blast_format_task(db, db_path, 'prot')
            )
            user_databases[os.path.basename(db)] = db_path

    tasks.append(
            get_link_file_task(os.path.abspath(transcriptome))
    )

    '''
    Calculate assembly information. First it runs some basic stats like N50 and
    number of contigs, and uses the HyperLogLog counter from khmer to
    estimate unique k-mers for checking redundancy. Then it runs BUSCO to
    assess completeness. These tasks are grouped under the 'assess' task.
    '''
    assess_tasks = []
    assess_tasks.append(
        get_transcriptome_stats_task(transcriptome, 
                                     os.path.basename(transcriptome + '.stats'))
    )
    results['stats'] = transcriptome + '.stats'

    '''
    BUSCO assesses completeness using a series of curated databases of core
    conserved genes.
    '''
    busco_cfg = common.CONFIG['settings']['busco']
    busco_output_name = '{0}.busco.results'.format(transcriptome)
    busco_dir = prog_paths.get('busco', '')
    assess_tasks.append(
        get_busco_task(transcriptome, busco_output_name, database_dict['BUSCO'],
                       'trans', n_threads, busco_cfg,
                       busco_dir=busco_dir)
    )
    results['busco'] = os.path.join('run_' + busco_output_name, 
                                    'short_summary_' + busco_output_name)

    # Collect the stats and BUSCO tasks under an "assess" group for convenience
    tasks.extend(assess_tasks)
    tasks.append(get_group_task('assess', assess_tasks))

    '''
    Run TransDecoder. TransDecoder first finds long ORFs with
    TransDecoder.LongOrfs, which are output as a FASTA file of protein
    sequences. We can then use these sequences to search against Pfam-A for
    conserved domains. TransDecoder.Predict uses the Pfam results to train its
    model for prediction of gene features.
    '''

    annotate_tasks = []

    transdecoder_output_dir = transcriptome + '.transdecoder_dir'
    transdecoder_dir = prog_paths.get('transdecoder', '')
    orf_cfg = common.CONFIG['settings']['transdecoder']['longorfs']
    annotate_tasks.append(
        get_transdecoder_orf_task(transcriptome, 
                                  orf_cfg,
                                  transdecoder_dir=transdecoder_dir)
    )
    orf_results = os.path.join(transdecoder_output_dir,
                               'longest_orfs.pep')
    results['ORFs'] = orf_results

    hmmer_dir = prog_paths.get('hmmer', '')
    pfam_results = transcriptome + '.pfam-A.tbl'
    annotate_tasks.append(
        get_hmmscan_task(orf_results, pfam_results,
                     database_dict['PFAM'], n_threads, 
                     common.CONFIG['settings']['hmmer']['hmmscan'],
                     hmmer_dir=hmmer_dir)
    )
    results['pfam'] = pfam_results

    predict_cfg = common.CONFIG['settings']['transdecoder']['predict']
    annotate_tasks.append(
        get_transdecoder_predict_task(transcriptome, 
                                      pfam_results,
                                      n_threads,
                                      predict_cfg,
                                      transdecoder_dir=transdecoder_dir)
    )
    protein_prediction_results = transcriptome + '.transdecoder.pep'
    results['prot_predictions'] = protein_prediction_results


    cmscan_cfg = common.CONFIG['settings']['infernal']['cmscan']
    infernal_dir = prog_paths.get('infernal', '')
    rfam_results = transcriptome + '.rfam.tbl'
    annotate_tasks.append(
        get_cmscan_task(transcriptome, rfam_results,
                     database_dict['RFAM'], n_threads, 
                     cmscan_cfg, infernal_dir=infernal_dir)
    )
    results['rfam'] = rfam_results

    '''
    Build the BLAST+ database for the transcriptome, which we'll need
    for doing RBH searches
    '''
    blast_dir = prog_paths.get('blast','')
    blast_cfg = common.CONFIG['settings']['blast']
    annotate_tasks.append(
        get_blast_format_task(transcriptome, transcriptome + '.db', 
                               'nucl', blast_dir=blast_dir)
    )
    
    last_dir = prog_paths.get('last', '')
    lastal_cfg = common.CONFIG['settings']['last']['lastal']
    orthodb_fn = database_dict['ORTHODB']
    tr_x_orthodb_fn = '{0}.x.orthodb.maf'.format(transcriptome)
    annotate_tasks.append(
        get_lastal_task(transcriptome, orthodb_fn, tr_x_orthodb_fn, True,
                       n_threads, lastal_cfg, last_dir=last_dir)
    )
    results['orthodb'] = tr_x_orthodb_fn
 
    results['user'] = {}
    for key, path in user_database_dict.iteritems():
        tr_x_db = '{0}.x.{1}.tsv'.format(transcriptome, key)
        db_x_tr = '{0}.x.{1}.tsv'.format(key, transcriptome)
        annotate_tasks.append(
            get_blast_task(transcriptome, path, 'blastx', tr_x_db,
                           n_threads, blast_cfg, blast_dir=blast_dir)
        )
        annotate_tasks.append(
            get_blast_task(path, transcriptome, 'tblastn', db_x_tr,
                           n_threads, blast_cfg, blast_dir=blast_dir)
        )
        results['user'][key] = (tr_x_db, db_x_tr)



    tasks.extend(annotate_tasks)

    return results, tasks

def run_annotate_tasks(transcriptome, output_dir, tasks, args=['run']):
    '''
    Set up doit's config for the actual analysis tasks.
    We'll put the doit database for these tasks into the output
    directory so that we don't end up scattering them around the
    filesystem, or worse, with one master db containing dependency
    metadata from every analysis ever run by the user!
    '''

    doit_config = {
                    'backend': common.DOIT_BACKEND,
                    'verbosity': common.DOIT_VERBOSITY,
                    'dep_file': os.path.join(output_dir, '.' +
                                             os.path.basename(transcriptome) +
                                             '.doit.db')
                  }

    cwd = os.getcwd()
    try:
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        os.chdir(output_dir)

        common.run_tasks(tasks, args, config=doit_config)
    finally:
        os.chdir(cwd)

