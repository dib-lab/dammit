#!/usr/bin/enb python
from __future__ import print_function

import logging
import os
from os import path
import sys

from doit.dependency import Dependency, SqliteDB

from . import common
from .log import LogReporter
from .hmmer import hmmpress
from .infernal import get_cmpress_task
from .last import get_lastdb_task
from .tasks import get_download_and_gunzip_task, \
                   get_download_and_untar_task, \
                   print_tasks

def get_handler(args, config):

    logger = logging.getLogger('DatabaseHandler')

    if args.database_dir is None:
        try:
            directory = os.environ['DAMMIT_DB_DIR']
            logger.debug('found DAMMIT_DB_DIR env variable')
        except KeyError:
            logger.debug('no DAMMIT_DB_DIR or --database-dir, using'\
                              ' default')
            directory = path.join(config['settings']['dammit_dir'],
                                  config['settings']['db_dir'])
    else:
        directory = args.database_dir
    directory = path.abspath(directory)
    doit_db = path.join(directory, 'databases.doit.db')

    handler = TaskHandler(doit_db, logger, config=config,
                          reporter=LogReporter(logger), 
                          backend=config['settings']['doit_backend'],
                          verbosity=args.verbosity,
                          continue=True)
    return handler


def register_builtin_tasks(handler, config, databases):

    settings = config['settings']

    register_pfam_tasks(handler, settings, databases)
    register_rfam_tasks(handler, settings, databases)
    register_orthodb_tasks(handler, settings, databases)
    register_busco_tasks(handler, settings, databases)
    register_uniref90_tasks(handler, settings, databases)
                          

def register_pfam_tasks(handler, settings, databases):
    pfam_A = databases['Pfam-A']
    handler.register_task('download:Pfam-A',
                          get_download_and_gunzip_task(pfam_A['url'],
                                                       pfam_A['filename']),
                          files={'Pfam-A': pfam_A['filename']})
    handler.register_task('hmmpress:Pfam-A',
                          get_hmmpress_task(pfam_A['filename'], settings['hmmer']))


def register_rfam_tasks(handler, settings, databases):
    rfam = databases['Rfam']
    handler.register_task('download:Rfam',
                           get_download_and_gunzip_task(rfam['url'], 
                                                        rfam['filename']),
                           files={'Rfam': rfam['filename']})
    handler.register_task('cmpress:Rfam',
                          get_cmpress_task(rfam['filename'], settings['infernal']))


def register_orthodb_tasks(handler, settings, databases):
    orthodb = databases['OrthoDB']
    handler.register_task('download:OrthoDB',
                          get_download_and_gunzip_task(orthodb['url'], 
                                                       orthodb['filename']),
                          files={'OrthoDB': orthodb['filename']})
    handler.register_task('lastdb:OrthoDB',
                          get_lastdb_task(orthodb['filename'], 
                                          orthodb['filename'], 
                                          settings['last']['lastdb'], 
                                          prot=True))


def register_busco_tasks(handler, settings, databases):
    busco = databases['busco']
    busco_dir = 'buscodb'
    for group_name in busco:
        group = busco[group_name]
        files = {'BUSCO-{0}'.format(group_name): path.join(busco_dir, group_name)}
        handler.register_task('download:BUSCO-{0}'.format(group_name),
                              get_downlod_and_untar_task(group['url'],
                                                         busco_dir,
                                                         label=group_name),
                              files=files)


def register_uniref90_tasks(handler, settings, databases):
    uniref90 = databases['uniref90']
    handler.register_task('download:uniref90',
                          get_download_and_gunzip_task(uniref90['url'],
                                                       uniref90['filename']),
                          files={'uniref90': uniref90['filename']})
    handler.register_task('lastdb:uniref90',
                          get_lastdb_task(uniref['filename'],
                                          uniref['filename'],
                                          settings['last']['lastdb'],
                                          prot=True))


def install(handler):
    uptodate, missing = handler.check_uptodate()
    if not uptodate:




class DatabaseHandler(object):

    def handle(self, doit_args=['run']):

        missing = self.check()
        
        if self.args.install:
            if missing:
                common.print_header('Installing databases', level=2)
                common.run_tasks(self.tasks, doit_args, config=self.doit_config)
            else:
                common.print_header('Nothing to install', level=2)
        else:
            if missing:
                sys.exit(1)

    def check_or_fail(self):
        missing = self.check()
        if missing:
            self.logger.error('Install databases to continue.')
            common.print_header('to prepare databases, run: dammit databases'\
                                ' --install', level=2)
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
            self.logger.warning('Database prep incomplete...')
        else:
            self.logger.info('All databases prepared!')

        return missing

