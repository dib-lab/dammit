#!/usr/bin/env python
from __future__ import print_function

import argparse
import logging
import os
import sys

from dammit import __version__
from dammit import common
from dammit import annotate
from dammit import databases
from dammit import dependencies
from dammit import log
from dammit.tasks import print_tasks


class DammitApp(object):

    def __init__(self, arg_src=sys.argv[1:]):
        print(arg_src)
        self.logger = logging.getLogger(self.__class__.__name__)
        self.meta = '{0}\n{1} {2}'.format(common.CONFIG['meta']['description'],
                                          ', '.join(common.CONFIG['meta']['authors']),
                                          common.CONFIG['meta']['date'])
        
        self.parser = self.get_parser()
        self.args = self.parser.parse_args(arg_src)

    def run(self):
        common.print_header(self.meta, 0)
        self.args.func()

    def get_parser(self):
        '''
        Build the main parser.
        '''
        parser = argparse.ArgumentParser(
                     description=common.add_header(self.meta, 0),
                     formatter_class=argparse.RawDescriptionHelpFormatter
                     )

        parser.add_argument('--debug', action='store_true', default=False)
        parser.add_argument('--version', action='version', 
                            version='%(prog)s ' +__version__)
        subparsers = parser.add_subparsers(title='dammit subcommands')

        def add_common_args(parser):
            ''' Add shared options to a parser.

            Shared options are added this way instead of to the main parser
            because I'd rather they come after the subcommand name.

            Args:
                parser (object): The parser to which arguments will be added.
            '''
            parser.add_argument('--database-dir', 
                                default=None, 
                                help='Directory to store databases. Existing'\
                                     ' databases will not be overwritten.'\
                                     ' By default, the database directory is'\
                                     ' $HOME/.dammit/databases.'
                                )

            parser.add_argument('--full', 
                                action='store_true', 
                                default=False,
                                help='Do complete run with uniref90. By default'\
                                     ' uniref90 is left out, as it is huge '\
                                     ' and homology searches take a long'\
                                     ' time.'
                                )

            parser.add_argument('--busco-group', 
                                default='metazoa',
                                choices=['metazoa', 'eukaryota', 'vertebrata',
                                         'arthropoda'],
                                help='Which BUSCO group to use. Depends on'\
                                     ' the organism being annotated.'
                                )

        '''
        Add the databases subcommand.
        '''
        desc = '''Check for databases and optionally download and prepare them 
               for use. By default, only check their status.'''
        databases_parser = subparsers.add_parser(
                               'databases',
                                description=desc,
                                help=desc
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
        Add the dependencies subcommand.
        '''
        desc = '''Checks for dependencies on system PATH. Unlike with the 
               databases, dependencies are not downloaded when missing and must 
               be installed by the user.'''
        dependencies_parser = subparsers.add_parser(
                                  'dependencies',
                                  description=desc,
                                  help=desc
                                  )
        dependencies_parser.set_defaults(func=self.handle_dependencies)

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
                              help=desc
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

        annotate_parser.add_argument('--n_threads', 
                                     type=int, 
                                     default=1,
                                     help='Number of threads to pass to programs'\
                                          ' supporting multithreading'
                                     )

        annotate_parser.add_argument('--user-databases', 
                                     nargs='+', 
                                     default=[],
                                     help='Optional additional protein databases. '\
                                          ' These will be searched with CRB-blast.'
                                     )

        add_common_args(annotate_parser)
        annotate_parser.set_defaults(func=self.handle_annotate)

        return parser

    def handle_databases(self):
        common.print_header('submodule: databases', level=1)

        dependencies.DependencyHandler().check_or_fail()

        databases.DatabaseHandler(self.args).handle()

        
    def handle_dependencies(self):
        common.print_header('submodule: dependencies', level=1)

        dependencies.DependencyHandler().handle()

    def handle_annotate(self):

        common.print_header('submodule: annotate', level=1)

        dependencies.DependencyHandler().check_or_fail()

        db_handler = databases.DatabaseHandler(self.args)
        db_handler.check_or_fail()

        annotate.AnnotateHandler(self.args, db_handler.databases).handle()


