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
                   get_remap_hmmer_task, \
                   get_cmscan_task, \
                   get_lastal_task, \
                   get_crb_blast_task, \
                   get_sanitize_fasta_task, \
                   get_rename_transcriptome_task, \
                   get_transeq_task, \
                   print_tasks

logger = logging.getLogger(__name__)

class AnnotateHandler(object):

    def __init__(self, args, database_dict):
        self.args = args
        self.logger = logging.getLogger(self.__class__.__name__)

        self._init_filenames()

        self.doit_config = {
                        'reporter': LogReporter(logger),
                        'backend': common.DOIT_BACKEND,
                        'verbosity': common.DOIT_VERBOSITY,
                        'continue': True,
                        'dep_file': os.path.join(self.directory, '.' +
                                                 os.path.basename(self.transcriptome_fn) +
                                                 '.doit.db')
                      }

        self.database_dict = database_dict
        self.tasks = list(self.get_tasks())

    def _init_filenames(self):
        '''Initialize all the input/output filename atributes.
        '''

        self.input_transcriptome_fn = os.path.abspath(self.args.transcriptome)
        self.transcriptome_fn = os.path.basename(self.input_transcriptome_fn)
        self.names_fn = '{0}.dammit.names.csv'.format(self.transcriptome_fn)
        if self.args.output_dir is None:
            out_dir = os.path.basename(self.input_transcriptome_fn) + '.dammit'
        else:
            out_dir = self.args.output_dir
        self.directory = os.path.abspath(out_dir)
        
        self.stats_fn = self.transcriptome_fn + '.stats.json'

        self.busco_basename = '{0}.{1}.busco.results'.format(self.transcriptome_fn,
                                                             self.args.busco_group)
        self.busco_dir = 'run_{0}'.format(self.busco_basename)
        busco_summary_fn = 'short_summary_{0}'.format(self.transcriptome_fn)
        self.busco_summary_fn = os.path.join(self.busco_dir,
                                             busco_summary_fn)
        
        self.translated_fn = '{0}.pep'.format(self.transcriptome_fn)

        self.transdecoder_dir = '{0}.transdecoder_dir'.format(self.transcriptome_fn)
        self.transdecoder_orf_fn = os.path.join(self.transdecoder_dir,
                                                'longest_orfs.pep')
        self.transdecoder_orf_gff3_fn = os.path.join(self.transdecoder_dir,
                                                     'longest_orfs.gff3')
        self.transdecoder_pfam_fn = '{0}.pfam.tbl'.format(self.transdecoder_orf_fn)
        self.transdecoder_pep_fn = '{0}.transdecoder.pep'.format(self.transcriptome_fn)
        self.transdecoder_gff3_fn = '{0}.transdecoder.gff3'.format(self.transcriptome_fn)
        
        self.pfam_fn = '{0}.pfam.csv'.format(self.transcriptome_fn)
        self.rfam_fn = '{0}.rfam.tbl'.format(self.transcriptome_fn)

        self.orthodb_fn = '{0}.x.orthodb.maf'.format(self.transcriptome_fn)

        self.user_pep_fn_dict = {}


    def handle(self):

        common.print_header('Running annotate!', level=2)
        self.logger.info('Transcriptome file: {0}'.format(self.transcriptome_fn))
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

            common.run_tasks(self.tasks, 
                             doit_args, 
                             config=self.doit_config)
        finally:
            self.logger.debug('chdir: {0}'.format(cwd))
            os.chdir(cwd)

    def rename_task(self):
        return get_rename_transcriptome_task(self.input_transcriptome_fn,
                                             self.transcriptome_fn,
                                             self.names_fn,
                                             self.args.name)

    def stats_task(self):
        '''Calculate assembly information. First it runs some basic stats like N50 and
        number of contigs, and uses the HyperLogLog counter from khmer to
        estimate unique k-mers for checking redundancy. Then it runs BUSCO to
        assess completeness. These tasks are grouped under the 'assess' task.
        '''

        return get_transcriptome_stats_task(self.transcriptome_fn, 
                                            self.stats_fn)
   
    def busco_task(self):
        '''BUSCO assesses completeness using a series of curated databases of core
        conserved genes.
        '''

        busco_cfg = common.CONFIG['settings']['busco']
        return get_busco_task(self.transcriptome_fn, 
                              self.busco_basename, 
                              self.database_dict['BUSCO'],
                              'trans', 
                              self.args.n_threads, 
                              busco_cfg)

    def transeq_task(self):
        '''Run transeq to do full six-frame protein translation.
        '''

        return get_transeq_task(self.transcriptome_fn,
                                self.translated_fn)

    def transdecoder_tasks(self):
        '''Run TransDecoder. TransDecoder first finds long ORFs with
        TransDecoder.LongOrfs, which are output as a FASTA file of protein
        sequences. We can then use these sequences to search against Pfam-A for
        conserved domains. TransDecoder.Predict uses the Pfam results to train its
        model for prediction of gene features.
        '''


        orf_cfg = common.CONFIG['settings']['transdecoder']['longorfs']
        yield get_transdecoder_orf_task(self.transcriptome_fn, 
                                        orf_cfg)

        yield get_hmmscan_task(self.transdecoder_orf_fn, 
                               self.transdecoder_pfam_fn,
                               self.database_dict['PFAM'], 
                               self.args.evalue,
                               self.args.n_threads, 
                               common.CONFIG['settings']['hmmer']['hmmscan'])

        yield get_remap_hmmer_task(self.transdecoder_pfam_fn,
                                   self.transdecoder_orf_gff3_fn,
                                   self.pfam_fn)

        predict_cfg = common.CONFIG['settings']['transdecoder']['predict']
        yield get_transdecoder_predict_task(self.transcriptome_fn, 
                                            self.transdecoder_pfam_fn,
                                            predict_cfg)

    def cmscan_task(self):
        '''Run Infernal. Infernal uses covariance models to detect
        RNA secondary structures. Here we use Rfam as our reference
        database.
        '''

        cmscan_cfg = common.CONFIG['settings']['infernal']['cmscan']
        return get_cmscan_task(self.transcriptome_fn, 
                               self.rfam_fn,
                               self.database_dict['RFAM'], 
                               self.args.evalue,
                               self.args.n_threads, 
                               cmscan_cfg)

    def orthodb_task(self):
        '''Run LAST to get homologies with OrthoDB. We use LAST here because
        it is much faster than BLAST+, and OrthoDB is pretty huge.
        '''
        
        lastal_cfg = common.CONFIG['settings']['last']['lastal']
        orthodb = self.database_dict['ORTHODB']
        return get_lastal_task(self.transcriptome_fn, 
                               orthodb, 
                               self.orthodb_fn, 
                               True,
                               self.args.evalue,
                               self.args.n_threads, 
                               lastal_cfg)

    def user_crb_tasks(self):
        '''Run conditional recipricol best hits LAST (CRBL) against the
        user-supplied databases.
        '''

        crb_blast_cfg = common.CONFIG['settings']['crb-blast']
        for path in self.args.user_databases:
            key = os.path.basename(path)
            yield get_sanitize_fasta_task(os.path.abspath(path),
                                          key)

            fn = '{0}.x.{1}.crbb.tsv'.format(self.transcriptome_fn, key)
            self.user_pep_fn_dict[key] = fn
            yield get_crb_blast_task(self.transcriptome_fn, 
                                     key, 
                                     fn, 
                                     self.args.evalue,
                                     crb_blast_cfg, 
                                     self.args.n_threads)

    def get_tasks(self):

        yield self.rename_task()
        yield self.stats_task()
        yield self.busco_task()
        yield self.transeq_task()
        for task in self.transdecoder_tasks():
            yield task
        yield self.cmscan_task()
        yield self.orthodb_task()
        for task in self.user_crb_tasks():
            yield task

        self.outputs, report_tasks = get_report_tasks(self.transcriptome_fn, 
                                                      self,
                                                      self.database_dict,
                                                      n_threads=self.args.n_threads)
        for task in report_tasks:
            yield task

