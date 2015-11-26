#!/usr/bin/env python
from __future__ import print_function

import logging
import os
from platform import system
import sys

from . import common
from .log import LogReporter
#from .crbl import CRBL
from .report import get_report_tasks
from .tasks import get_transcriptome_stats_task, \
                   get_busco_task, \
                   get_group_task, \
                   get_link_file_task, \
                   get_transdecoder_predict_task, \
                   get_transdecoder_orf_task, \
                   get_hmmscan_task, \
                   get_cmscan_task, \
                   get_lastal_task, \
                   get_crb_blast_task, \
                   get_sanitize_fasta_task, \
                   print_tasks

logger = logging.getLogger(__name__)

class AnnotateHandler(object):

    def __init__(self, args, database_dict):
        self.input_transcriptome = os.path.abspath(args.transcriptome)
        self.transcriptome = os.path.basename(self.input_transcriptome)
        self.args = args
        self.logger = logging.getLogger(self.__class__.__name__)

        if args.output_dir is None:
            out_dir = os.path.basename(self.input_transcriptome) + '.dammit'
        else:
            out_dir = args.output_dir
        self.directory = os.path.abspath(out_dir)


        self.doit_config = {
                        'reporter': LogReporter(logger),
                        'backend': common.DOIT_BACKEND,
                        'verbosity': common.DOIT_VERBOSITY,
                        'continue': True,
                        'dep_file': os.path.join(self.directory, '.' +
                                                 os.path.basename(self.transcriptome) +
                                                 '.doit.db')
                      }

        self.database_dict = database_dict
        self.results, self.tasks = self.get_tasks()

    def handle(self):

        common.print_header('Running annotate!', level=2)
        self.logger.info('Transcriptome file: {0}'.format(self.transcriptome))
        self.logger.info('Output directory: {0}'.format(self.directory))
        print_tasks(self.tasks, logger=self.logger)

        self.run_tasks()

    def run_tasks(self, doit_args=['run']):
        '''
        Set up doit's config for the actual analysis tasks.
        We'll put the doit database for these tasks into the output
        directory so that we don't end up scattering them around the
        filesystem, or worse, with one master db containing dependency
        metadata from every analysis ever run by the user!
        '''

        cwd = os.getcwd()
        self.logger.debug('cwd: {0}'.format(cwd))
        try:
            if not os.path.exists(self.directory):
                self.logger.debug('makedirs: {0}'.format(self.directory))
                os.makedirs(self.directory)
            os.chdir(self.directory)

            common.run_tasks(self.tasks, doit_args, config=self.doit_config)
        finally:
            self.logger.debug('chdir: {0}'.format(cwd))
            os.chdir(cwd)


    def get_tasks(self):

        tasks = []
        results = {}

        tasks.append(
                get_sanitize_fasta_task(self.input_transcriptome,
                                        self.transcriptome)
        )

        '''
        Calculate assembly information. First it runs some basic stats like N50 and
        number of contigs, and uses the HyperLogLog counter from khmer to
        estimate unique k-mers for checking redundancy. Then it runs BUSCO to
        assess completeness. These tasks are grouped under the 'assess' task.
        '''
        assess_tasks = []
        assess_tasks.append(
            get_transcriptome_stats_task(self.transcriptome, 
                                         os.path.basename(self.transcriptome + '.stats'))
        )
        results['stats'] = self.transcriptome + '.stats'

        '''
        BUSCO assesses completeness using a series of curated databases of core
        conserved genes.
        '''
        busco_cfg = common.CONFIG['settings']['busco']
        busco_output_name = '{0}.busco.results'.format(self.transcriptome)
        assess_tasks.append(
            get_busco_task(self.transcriptome, busco_output_name, self.database_dict['BUSCO'],
                           'trans', self.args.n_threads, busco_cfg)
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

        transdecoder_output_dir = self.transcriptome + '.transdecoder_dir'
        orf_cfg = common.CONFIG['settings']['transdecoder']['longorfs']
        annotate_tasks.append(
            get_transdecoder_orf_task(self.transcriptome, 
                                      orf_cfg)
        )
        orf_pep = os.path.join(transdecoder_output_dir,
                                   'longest_orfs.pep')
        orf_gff3 = os.path.join(transdecoder_output_dir,
                                'longest_orfs.gff3')
        results['ORFs_pep'] = orf_pep
        results['ORFs_gff3'] = orf_gff3

        pfam_results = self.transcriptome + '.pfam-A.tbl'
        annotate_tasks.append(
            get_hmmscan_task(orf_pep, pfam_results,
                         self.database_dict['PFAM'], self.args.n_threads, 
                         common.CONFIG['settings']['hmmer']['hmmscan'])
        )
        results['pfam'] = pfam_results

        predict_cfg = common.CONFIG['settings']['transdecoder']['predict']
        annotate_tasks.append(
            get_transdecoder_predict_task(self.transcriptome, 
                                          pfam_results,
                                          predict_cfg)
        )
        protein_prediction_pep = self.transcriptome + '.transdecoder.pep'
        protein_prediction_gff3 = self.transcriptome + '.transdecoder.gff3'
        results['prot_predictions_pep'] = protein_prediction_pep
        results['prot_predictions_gff3'] = protein_prediction_gff3

        '''
        Run Infernal. Infernal uses covariance models to detect
        RNA secondary structures. Here we use Rfam as our reference
        database.
        '''
        cmscan_cfg = common.CONFIG['settings']['infernal']['cmscan']
        rfam_results = self.transcriptome + '.rfam.tbl'
        annotate_tasks.append(
            get_cmscan_task(self.transcriptome, rfam_results,
                         self.database_dict['RFAM'], self.args.n_threads, 
                         cmscan_cfg)
        )
        results['rfam'] = rfam_results

        '''
        Run LAST to get homologies with OrthoDB. We use LAST here because
        it is much faster than BLAST+, and OrthoDB is pretty huge.
        '''
        
        lastal_cfg = common.CONFIG['settings']['last']['lastal']
        orthodb_fn = self.database_dict['ORTHODB']
        tr_x_orthodb_fn = '{0}.x.orthodb.maf'.format(self.transcriptome)
        annotate_tasks.append(
            get_lastal_task(self.transcriptome, orthodb_fn, tr_x_orthodb_fn, True,
                           self.args.n_threads, lastal_cfg)
        )
        results['orthodb'] = tr_x_orthodb_fn
     

        '''
        Run conditional recipricol best hits LAST (CRBL) against the
        user-supplied databases.
        '''
        results['user'] = {}
        crb_blast_cfg = common.CONFIG['settings']['crb-blast']
        for path in self.args.user_databases:
            key = os.path.basename(path)
            tasks.append(
                get_sanitize_fasta_task(os.path.abspath(path),
                                        key)
            )
            fn = '{0}.x.{1}.crbb.tsv'.format(self.transcriptome, key)
            annotate_tasks.append(
                get_crb_blast_task(self.transcriptome, key, fn, 
                                   crb_blast_cfg, self.args.n_threads)
            )
            results['user'][key] = fn

        tasks.extend(annotate_tasks)

        outputs, report_tasks = get_report_tasks(self.transcriptome, results,
                                                 self.database_dict,
                                                 n_threads=self.args.n_threads)
        tasks.extend(report_tasks)

        return results, tasks

