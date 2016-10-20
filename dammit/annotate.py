#!/usr/bin/env python
from __future__ import print_function

from collections import OrderedDict
import logging
import os
from os import path
import sys

from shmlast.app import CRBL
from shmlast.last import lastal_task as get_lastal_task

from .handler import TaskHandler

from .tasks.fastx import (get_transcriptome_stats_task,
                                get_sanitize_fasta_task,
                                get_rename_transcriptome_task,
                                get_annotate_fasta_task)
from .tasks.busco import get_busco_task
from .tasks.utils import get_group_task
from .tasks.shell import get_link_file_task
from .tasks.gff import (get_maf_gff3_task,
                        get_shmlast_gff3_task,
                        get_hmmscan_gff3_task,
                        get_cmscan_gff3_task,
                        get_gff3_merge_task,
                        get_maf_best_hits_task)
from .tasks.hmmer import get_hmmscan_task, get_remap_hmmer_task
from .tasks.infernal import get_cmscan_task
from .tasks.transdecoder import (get_transdecoder_predict_task,
                                       get_transdecoder_orf_task)
from . import ui
from . import log

logger = logging.getLogger(__name__)

def get_handler(config, databases, builtins=True):

    logger = logging.getLogger('AnnotateHandler')

    if config['output_dir'] is None:
        out_dir = path.basename(config['transcriptome'] + '.dammit')
    else:
        out_dir = config['output_dir']
    directory = path.abspath(out_dir)
    
    handler = TaskHandler(directory, logger, config=config,
                          db='annotate',
                          backend=config['doit_backend'],
                          verbosity=config['verbosity'])
    log.start_logging(path.join(directory, 'dammit.log'))

    input_fn = path.join(directory, path.basename(config['transcriptome']))
    name_map_fn = input_fn + '.dammit.namemap.csv'
    handler.register_task('rename-transcriptome',
                          get_rename_transcriptome_task(path.abspath(config['transcriptome']),
                                                        input_fn,
                                                        name_map_fn,
                                                        config['name']),
                          files={'transcriptome': input_fn,
                                 'name_map': name_map_fn})
    
    if builtins:
        return register_builtin_tasks(handler, config, databases)
    return handler


def run_annotation(handler):
    print(ui.header('Annotation', level=3))
    print('Doit Database: {0}'.format(handler.dep_file))
    print('Input Transcriptome: {0}'.format(handler.files['transcriptome']))
    msg = '*All annotation tasks up-to-date.*'
    uptodate, statuses = handler.print_statuses(uptodate_msg=msg)
    if not uptodate:
        print('Running pipeline...')
        handler.run()
    else:
        print('Pipeline is already completed!')
        sys.exit(0)


def register_builtin_tasks(handler, config, databases):

    register_stats_task(handler)
    register_busco_task(handler, config, databases)
    register_transdecoder_tasks(handler, config, databases)
    register_rfam_tasks(handler, config, databases)
    register_last_tasks(handler, config, databases)
    register_user_db_tasks(handler, config, databases)
    register_annotate_tasks(handler, config, databases)

    return handler


def register_stats_task(handler):
    input_fn = handler.files['transcriptome']
    stats_fn = input_fn + '.dammit.stats.json'
    handler.register_task('transcriptome-stats',
                          get_transcriptome_stats_task(input_fn,
                                                       stats_fn),
                          files={'stats': stats_fn})

def register_busco_task(handler, config, databases):
    input_fn = handler.files['transcriptome']
    busco_group = config['busco_group']
    busco_database = databases['BUSCO-{0}'.format(busco_group)]
    busco_basename = '{0}.{1}.busco.results'.format(input_fn, busco_group)
    busco_out_dir = 'run_{0}'.format(busco_basename)

    handler.register_task('BUSCO-{0}'.format(busco_group),
                          get_busco_task(input_fn,
                                         busco_basename,
                                         busco_database,
                                         input_type='trans',
                                         n_threads=config['n_threads'],
                                         params=config['busco']['params']),
                          files={'BUSCO': busco_out_dir})


def register_transdecoder_tasks(handler, config, databases):
    '''Run TransDecoder. TransDecoder first finds long ORFs with
    TransDecoder.LongOrfs, which are output as a FASTA file of protein
    sequences. We can then use these sequences to search against Pfam-A for
    conserved domains. TransDecoder.Predict uses the Pfam results to train its
    model for prediction of gene features.
    '''
    input_fn = handler.files['transcriptome']
    transdecoder_dir = '{0}.transdecoder_dir'.format(input_fn)
    
    handler.register_task('TransDecoder.LongOrfs',
                          get_transdecoder_orf_task(input_fn, 
                                                    params=config['transdecoder']['longorfs']),
                          files={'longest_orfs': path.join(transdecoder_dir, 'longest_orfs.pep')})
    pfam_fn = path.join(transdecoder_dir, 'longest_orfs.pep.x.pfam.tbl')
    handler.register_task('hmmscan:Pfam-A',
                          get_hmmscan_task(handler.files['longest_orfs'],
                                           pfam_fn,
                                           databases['Pfam-A'],
                                           cutoff=config['evalue'],
                                           n_threads=config['n_threads'],
                                           pbs=config['pbs'],
                                           params=config['hmmer']['hmmscan']),
                          files={'longest_orfs_pfam': pfam_fn})

    pfam_csv_fn = '{0}.x.pfam-A.csv'.format(input_fn)
    handler.register_task('hmmscan:Pfam-A:remap',
                          get_remap_hmmer_task(handler.files['longest_orfs_pfam'],
                                               path.join(transdecoder_dir, 'longest_orfs.gff3'),
                                               pfam_csv_fn),
                          files={'Pfam-A-csv': pfam_csv_fn})

    predict_cfg = config['transdecoder']['predict']
    transdecoder_pep = '{0}.transdecoder.pep'.format(input_fn)
    transdecoder_gff3 = '{0}.transdecoder.gff3'.format(input_fn)
    handler.register_task('TransDecoder.Predict',
                          get_transdecoder_predict_task(input_fn, 
                                                        pfam_fn,
                                                        predict_cfg),
                          files={'transdecoder-pep': transdecoder_pep,
                                 'transdecoder-gff3': transdecoder_gff3})

    pfam_gff3 = '{0}.x.pfam-A.gff3'.format(input_fn)
    handler.register_task('gff3:Pfam-A',
                           get_hmmscan_gff3_task(pfam_csv_fn,
                                                 pfam_gff3, 
                                                 'Pfam'),
                           files={'Pfam-A-gff3': pfam_gff3})


def register_rfam_tasks(handler, config, databases):
    input_fn = handler.files['transcriptome']
    output_fn = '{0}.x.rfam.tbl'.format(input_fn)
    handler.register_task('cmscan:Rfam',
                          get_cmscan_task(input_fn,
                                          output_fn,
                                          databases['Rfam'],
                                          cutoff=config['evalue'],
                                          n_threads=config['n_threads'],
                                          pbs=config['pbs'],
                                          params=config['infernal']['cmscan']),
                          files={'Rfam': output_fn})

    rfam_gff3 = '{0}.x.rfam.gff3'.format(input_fn)
    handler.register_task('gff3:Rfam',
                          get_cmscan_gff3_task(output_fn,
                                               rfam_gff3, 
                                               'Rfam'),
                          files={'Rfam-gff3': rfam_gff3})


def register_last_tasks(handler, config, databases):
    input_fn = handler.files['transcriptome']
    lastal_cfg = config['last']['lastal']
    
    dbs = OrderedDict()
    dbs['OrthoDB'] = databases['OrthoDB']
    if config['full'] is True:
        dbs['uniref90'] = databases['uniref90']

    for name, db in dbs.items():
        output_fn = '{0}.x.{1}.maf'.format(input_fn, name)
        handler.register_task('lastal:{0}'.format(name),
                              get_lastal_task(input_fn,
                                              db,
                                              output_fn,
                                              translate=True,
                                              cutoff=config['evalue'],
                                              n_threads=config['n_threads'],
                                              frameshift=lastal_cfg['frameshift'],
                                              pbs=config['pbs'],
                                              params=lastal_cfg['params']),
                              files={name: output_fn})

        best_fn = '{0}.x.{1}.best.csv'.format(input_fn, name)
        gff3_fn = '{0}.x.{1}.best.gff3'.format(input_fn, name)

        handler.register_task('lastal:best-hits:{0}'.format(name),
                              get_maf_best_hits_task(output_fn,
                                                     best_fn),
                              files={'{0}-best-hits'.format(name): best_fn})
        handler.register_task('gff3:{0}'.format(name),
                              get_maf_gff3_task(best_fn,
                                                gff3_fn,
                                                name),
                              files={'{0}-gff3'.format(name): gff3_fn})


def register_annotate_tasks(handler, config, databases):
    input_fn = handler.files['transcriptome']
    gff3_files = [fn for name, fn in handler.files.items() if name.endswith('-gff3')]
    merged_gff3 = '{0}.dammit.gff3'.format(input_fn)
    handler.register_task('gff3:merge-all',
                          get_gff3_merge_task(gff3_files, merged_gff3),
                          files={'merged-gff3': merged_gff3})

    annotated_fn = '{0}.dammit.fasta'.format(input_fn)
    handler.register_task('annotate:fasta',
                          get_annotate_fasta_task(input_fn,
                                                  merged_gff3,
                                                  annotated_fn),
                          files={'annotated-fasta': annotated_fn})


def register_user_db_tasks(handler, config, databases):
    '''Run conditional recipricol best hits LAST (CRBL) against the
    user-supplied databases.
    '''
    
    if not 'user_databases' in config:
        return

    input_fn = handler.files['transcriptome']
    for db_path in config['user_databases']:
        db_path = path.abspath(db_path)
        db_basename = path.basename(db_path)

        results_fn = '{0}.x.{1}.crbl.csv'.format(input_fn, db_basename)
        gff3_fn = '{0}.x.{1}.crbl.gff3'.format(input_fn, db_basename)

        crbl = CRBL(input_fn,
                    db_path, 
                    results_fn, 
                    n_threads=config['n_threads'],
                    cutoff=config['evalue'])
        for task in crbl.tasks():
            task.name = 'user-database-shmlast:{0}'.format(task.name)
            handler.register_task(task.name, task)
        handler.register_task('gff3:{0}'.format(results_fn),
                              get_shmlast_gff3_task(results_fn,
                                                    gff3_fn,
                                                    db_basename),
                              files={'{0}-crbl-gff3'.format(db_basename): gff3_fn})
        handler.files['{0}-crbl'.format(db_basename)] = results_fn

