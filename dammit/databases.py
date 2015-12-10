#!/usr/bin/enb python
from __future__ import print_function

import logging
import os
import sys

from doit.dependency import Dependency, SqliteDB

from . import common
from .log import LogReporter
from .tasks import get_download_and_gunzip_task, \
                   get_hmmpress_task, \
                   get_cmpress_task, \
                   get_download_and_untar_task, \
                   get_lastdb_task, \
                   print_tasks


class DatabaseHandler(object):

    def __init__(self, args):

        self.args = args
        self.logger = logging.getLogger(self.__class__.__name__)

        if args.database_dir is None:
            try:
                directory = os.environ['DAMMIT_DB_DIR']
                self.logger.debug('found DAMMIT_DB_DIR env variable')
            except KeyError:
                self.logger.debug('no DAMMIT_DB_DIR or --database-dir, using'\
                                  ' default')
                directory = os.path.join(common.get_dammit_dir(), 
                                         common.CONFIG['settings']['db_dir'])
        else:
            directory = args.database_dir
        self.directory = os.path.abspath(directory)

        self.doit_config = {
                            'reporter': LogReporter(self.logger),
                            'backend': common.DOIT_BACKEND,
                            'verbosity': common.DOIT_VERBOSITY,
                            'continue': True,
                            'dep_file': os.path.join(self.directory, 'databases.doit.db')
                           }
        self.logger.debug('doit_config:{0}'.format(self.doit_config))

        self.logger.debug('database dir: {0}'.format(self.directory))
        try:
            os.mkdir(self.directory)
        except OSError:
            self.logger.debug('database dir already exists')

        self.databases, self.tasks = self.get_tasks()


    def handle(self, doit_args=['run']):

        missing = self.check()
        print_tasks(self.tasks, logger=self.logger)
        
        if self.args.install:
            if missing:
                common.print_header('Installing databases', level=2)
                common.run_tasks(self.tasks, doit_args, config=self.doit_config)
            else:
                common.print_header('Nothing to install', level=2)

    def check_or_fail(self):
        missing = self.check()
        if missing:
            self.logger.error('Install databases to continue; exiting')
            sys.exit(1)

    def check(self):

        common.print_header('Checking for database prep (dir: {0})'.format(self.directory),
                            level=2)

        dep_manager = Dependency(SqliteDB, self.doit_config['dep_file'])
        missing = False
        for task in self.tasks:
            status = dep_manager.get_status(task, self.tasks)
            self.logger.debug('{0}:{1}'.format(task.name, status.status))
            if status.status != 'up-to-date':
                missing = True
                self.logger.warning('[ ] {0}'.format(task.name))
            else:
                self.logger.info('[x] {0}'.format(task.name))

        common.print_header('Database results', level=2)

        if missing:
            self.logger.warning('Database prep incomplete')
            common.print_header('to prepare databases, run: dammit databases'\
                                ' --install', level=2)
        else:
            self.logger.info('All databases prepared!')

        return missing

    def get_tasks(self):
        '''Generate tasks for installing the bundled databases. 
        
        These tasks download the databases, unpack them, and format them for use.
        Current bundled databases are:

            * Pfam-A (protein domans)
            * Rfam (RNA models)
            * OrthoDB8 (conserved ortholog groups)
            * uniref90 (protiens, if --full selected)
        
        User-supplied databases are downloaded separately.

        Args:
            self.directory (str): Directory where the databases will be stored.
            busco_db (str): The BUSCO group to use.
            full (bool): Whether to do a full run and get UNIREF90 as well.

        Returns:
            dict: A dictionary of the final database paths.
            list: A list of the doit tasks.

        '''

        tasks = []
        databases = {}

        # Get Pfam-A and prepare it for use with hmmer
        PFAM = os.path.join(self.directory, common.DATABASES['pfam']['filename'])
        tasks.append(
            get_download_and_gunzip_task(common.DATABASES['pfam']['url'], PFAM)
        )
        tasks.append(
            get_hmmpress_task(PFAM, common.CONFIG['settings']['hmmer'])
        )
        databases['PFAM'] = os.path.abspath(PFAM)

        # Get Rfam and prepare it for use with Infernal
        RFAM = os.path.join(self.directory, common.DATABASES['rfam']['filename'])
        tasks.append(
            get_download_and_gunzip_task(common.DATABASES['rfam']['url'], RFAM)
        )
        tasks.append(
            get_cmpress_task(RFAM, common.CONFIG['settings']['infernal'])
        )
        databases['RFAM'] = os.path.abspath(RFAM)

        # Get OrthoDB and prepare it for BLAST use
        ORTHODB = os.path.join(self.directory, common.DATABASES['orthodb']['filename'])
        tasks.append(
            get_download_and_gunzip_task(common.DATABASES['orthodb']['url'], ORTHODB)
        )

        lastdb_cfg = common.CONFIG['settings']['last']['lastdb']
        tasks.append(
            get_lastdb_task(ORTHODB, ORTHODB + '.db', lastdb_cfg, prot=True)
        )
        ORTHODB += '.db'
        databases['ORTHODB'] = os.path.abspath(ORTHODB)

        ORTHODB_GENES = os.path.join(self.directory,
                                     common.DATABASES['orthodb_genes']['filename'])
        tasks.append(
            get_download_and_gunzip_task(common.DATABASES['orthodb_genes']['url'],
                                         ORTHODB_GENES)
        )
        databases['ORTHODB_GENES'] = os.path.abspath(ORTHODB_GENES)

        # A little confusing. First, we get the top-level BUSCO path:
        BUSCO = os.path.join(self.directory, 'buscodb')
        tasks.append(
            # That top-level path is given to the download task:
            get_download_and_untar_task(common.DATABASES['busco'][self.args.busco_group]['url'], 
                                        BUSCO,
                                        label=self.args.busco_group)
        )
        # The untarred arhive has a folder named after the group:
        databases['BUSCO'] = os.path.abspath(os.path.join(BUSCO, self.args.busco_group))

        # Get uniref90 if the user specifies
        # Ignoring this until we have working CRBL
        if self.args.full and False:
            UNIREF = os.path.join(common.DATABASES['uniref90']['filename'])
            tasks.append(
                get_download_and_gunzip_task(common.DATABASES['uniref90']['url'], UNIREF)
            )
            tasks.append(
                get_blast_format_task(UNIREF, UNIREF + '.db',
                                      common.DATABASES['uniref90']['db_type'])
            )
            UNIREF += '.db'
            databases['UNIREF'] = os.path.abspath(UNIREF)

        return databases, tasks

