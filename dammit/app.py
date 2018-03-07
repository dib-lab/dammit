# Copyright (C) 2015-2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the BSD license.  See the LICENSE file for details.

import argparse
import glob
import logging
import os
import sys

from dammit import annotate
from dammit import databases
from dammit import log
from dammit import utils
from dammit import ui
from dammit.meta import __version__, __authors__, __description__, __date__, get_config
from dammit.annotate import (build_quick_pipeline,
                             build_default_pipeline,
                             build_full_pipeline)


class DammitApp(object):

    def __init__(self, arg_src=sys.argv[1:]):
        self.logger = logging.getLogger(self.__class__.__name__)

        self.config_d, self.databases_d = get_config()
        self.parser = self.get_parser()

        self.args = self.parser.parse_args(arg_src)
        if hasattr(self.args, 'config_file') and self.args.config_file is not None:
            with open(self.args.config_file) as fp:
                self.config_d.update(json.load(fp))
        self.config_d.update(vars(self.args))
 
    def run(self):
        print(ui.header('dammit'))
        print(ui.header(__description__, level=2))
        about = '\nby {0}\n\n**v{1}**, {2}\n'.format(', '.join(__authors__),
                                           __version__, __date__)
        print(about)
        return self.args.func()

    def description(self):
        return ui.header('dammit: ' + __description__)

    def epilog(self):
        return 'Available BUSCO groups are: '\
               '{0}'.format(', '.join(sorted(self.databases_d['BUSCO'].keys())))

    def get_parser(self):
        '''
        Build the main parser.
        '''
        parser = argparse.ArgumentParser(
                 description=self.description(),
                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.set_defaults(func=parser.print_help)

        parser.add_argument('--debug', action='store_true', default=False)
        parser.add_argument('--version', action='version',
                            version='%(prog)s ' + __version__)
        subparsers = parser.add_subparsers(title='dammit subcommands')

        def add_common_args(parser):
            ''' Add shared options to a parser.

            Shared options are added this way instead of to the main parser
            because I'd rather they come after the subcommand name.

            Args:
                parser (object): The parser to which arguments will be added.
            '''
            parser.add_argument('--database-dir',
                                default=databases.default_database_dir(self.logger),
                                help='Directory to store databases. Existing'\
                                     ' databases will not be overwritten.'\
                                     ' By default, the database directory is'\
                                     ' $HOME/.dammit/databases.'
                                )

            parser.add_argument('--busco-group',
                                default='metazoa',
                                metavar='[metazoa, eukaryota, vertebrata, ...]',
                                choices=list(self.databases_d['BUSCO'].keys()),
                                help='Which BUSCO group to use. Should be chosen'\
                                     ' based on the organism being annotated.'\
                                     ' Full list of options is below.'
                                )

            parser.add_argument('--n_threads',
                                type=int,
                                default=1,
                                help='For annotate, number of threads to pass to '\
                                     'programs  supporting multithreading. For '\
                                     'databases, number of simultaneous tasks '\
                                     'to execute.'
                                )

            parser.add_argument('--config-file',
                                help='A JSON file providing values to override'\
                                     ' built-in config. Advanced use only!'
                                )

            parser.add_argument('--verbosity',
                                default=0,
                                type=int,
                                choices=[0,1,2],
                                help='Verbosity level for doit tasks.'
                                )

            parser.add_argument('--profile',
                                default=False,
                                action='store_true',
                                help='Profile task execution.')
            
            parser.add_argument('--force',
                                default=False,
                                action='store_true',
                                help='Ignore missing database tasks.')

            pgroup = parser.add_mutually_exclusive_group()
            pgroup.add_argument('--full',
                                action='store_true',
                                default=False,
                                help='Run a "complete" annotation; includes'\
                                     ' uniref90, which is left out of the'\
                                     ' default pipeline because it is huge'\
                                     ' and homology searches take a long'\
                                     ' time.'
                                )

            pgroup.add_argument('--quick',
                                default=False,
                                action='store_true',
                                help='Run a "quick" annotation; excludes'\
                                     ' the Infernal Rfam tasks, the HMMER'\
                                     ' Pfam tasks, and the LAST OrthoDB'\
                                     ' and uniref90 tasks. Best for users'\
                                     ' just looking to get basic stats'\
                                     ' and conditional reciprocal best'\
                                     ' LAST from a protein database.')

        migrate_parser= subparsers.add_parser('migrate')
        migrate_parser.add_argument('--destructive', default=False,
                                    action='store_true')
        add_common_args(migrate_parser)
        migrate_parser.set_defaults(func=self.handle_migrate)


        '''
        Add the databases subcommand.
        '''
        desc = '''Check for databases and optionally download and prepare them
               for use. By default, only check their status.'''
        databases_parser = subparsers.add_parser(
                               'databases',
                                description=desc,
                                epilog=self.epilog(),
                                help=desc,
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter
                                )

        databases_parser.add_argument('--install',
                                      action='store_true',
                                      default=False,
                                      help='Install missing databases. Downloads'
                                           ' and preps where necessary'
                                      )

        add_common_args(databases_parser)
        databases_parser.set_defaults(func=self.handle_databases)

        '''
        Add the annotation subcommand.
        '''
        desc = '''The main annotation pipeline. Calculates assembly stats;
               runs BUSCO; runs LAST against OrthoDB (and optionally uniref90),
               HMMER against Pfam, Inferal against Rfam, and Conditional Reciprocal
               Best-hit Blast against user databases; and aggregates all results in
               a properly formatted GFF3 file.'''
        annotate_parser = subparsers.add_parser(
                              'annotate',
                              usage='%(prog)s <transcriptome> [OPTIONS]',
                              description=desc,
                              epilog=self.epilog(),
                              help=desc,
                              formatter_class=argparse.ArgumentDefaultsHelpFormatter
                              )

        annotate_parser.add_argument('transcriptome',
                                     help='FASTA file with the transcripts to be'\
                                          ' annotated.'
                                     )

        annotate_parser.add_argument('-n', '--name',
                                     default='Transcript',
                                     help='Base name to use for renaming the'\
                                          ' input transcripts. The new names'\
                                          ' will be of the form <name>_<X>.'\
                                          ' It should not have spaces, pipes,'\
                                          ' ampersands, or other characters'\
                                          ' with special meaning to BASH.'
                                     )

        annotate_parser.add_argument('-e', '--evalue',
                                     default=1e-5,
                                     type=float,
                                     help='e-value cutoff for similarity'\
                                          ' searches.'
                                     )

        annotate_parser.add_argument('-o', '--output-dir',
                                     default=None,
                                     help='Output directory. By default this will'\
                                          ' be the name of the transcriptome file'\
                                          ' with `.dammit` appended'
                                     )

        annotate_parser.add_argument('--user-databases',
                                     nargs='+',
                                     default=[],
                                     help='Optional additional protein databases. '\
                                          ' These will be searched with CRB-blast.'
                                     )

        annotate_parser.add_argument('--sshloginfile',
                                     default=None,
                                     help='Distribute execution across the specified nodes.')
                                    

        add_common_args(annotate_parser)
        annotate_parser.set_defaults(func=self.handle_annotate)

        return parser

    def handle_migrate(self):
        with utils.Move(self.args.database_dir):
            odb_files = glob.glob('aa_seq_euk.fasta.db.*')
            for fn in odb_files:
                pre, _, suf = fn.partition('.db')
                newfn = pre + suf
                if self.args.destructive:
                    os.rename(fn, newfn)
                else:
                    os.symlink(fn, newfn)

    def handle_databases(self):
        log.start_logging()
        print(ui.header('submodule: databases', level=2))

        handler = databases.get_handler(self.config_d)
        if self.args.quick:
            databases.build_quick_pipeline(handler,
                                           self.config_d,
                                           self.databases_d)
        else:
            databases.build_default_pipeline(handler, 
                                             self.config_d,
                                             self.databases_d,
                                             with_uniref=self.args.full)
        if self.args.install:
            return databases.install(handler)
        else:
            databases.check_or_fail(handler)

    def handle_annotate(self):
        log.start_logging()
        print(ui.header('submodule: annotate', level=2))

        db_handler = databases.get_handler(self.config_d)

        if self.args.quick:
            databases.build_quick_pipeline(db_handler,
                                           self.config_d,
                                           self.databases_d)
        else:
            databases.build_default_pipeline(db_handler, 
                                             self.config_d,
                                             self.databases_d,
                                             with_uniref=self.args.full)
        if self.config_d['force'] is True:
            utd_msg = '*All database tasks up-to-date.*'
            ood_msg = '*Some database tasks out-of-date; '\
                      'FORCE is True, ignoring!'
            uptodate, statuses = db_handler.print_statuses(uptodate_msg=utd_msg,
                                                           outofdate_msg=ood_msg)
        else:
            databases.check_or_fail(db_handler)

        annotate_handler = annotate.get_handler(self.config_d, db_handler.files)
        if self.args.quick:
            build_quick_pipeline(annotate_handler, 
                                 self.config_d, 
                                 db_handler.files)
        elif self.args.full:
            build_full_pipeline(annotate_handler, 
                                self.config_d, 
                                db_handler.files)
        else:
            build_default_pipeline(annotate_handler, 
                                   self.config_d, 
                                   db_handler.files)
        return annotate.run_annotation(annotate_handler)
