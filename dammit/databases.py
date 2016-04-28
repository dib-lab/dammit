#!/usr/bin/enb python
from __future__ import print_function

import logging
import os
from os import path
import sys

from doit.dependency import Dependency, SqliteDB

import .ui
from .handler import TaskHandler
from .log import LogReporter
from .hmmer import hmmpress
from .infernal import get_cmpress_task
from .last import get_lastdb_task
from .tasks import get_download_and_gunzip_task, \
                   get_download_and_untar_task, \
                   print_tasks

def get_handler(args, config, databases):

    logger = logging.getLogger('DatabaseHandler')

    if args.database_dir is None:
        try:
            directory = os.environ['DAMMIT_DB_DIR']
            logger.debug('found DAMMIT_DB_DIR env variable')
        except KeyError:
            logger.debug('no DAMMIT_DB_DIR or --database-dir, using'\
                              ' default')
            directory = path.join(os.environ['HOME'], '.dammit', 'databases')
    else:
        directory = args.database_dir
    directory = path.abspath(directory)
    doit_db = path.join(directory, 'databases.doit.db')

    handler = TaskHandler(doit_db, logger, config=config,
                          backend=config['settings']['doit_backend'],
                          verbosity=args.verbosity,
                          continue=True)

    return register_builtin_tasks(handler, config, databases)


def install(handler):
    uptodate, missing = handler.check_uptodate()
    if not uptodate:
        handler.print_outofdate()
        print('Installing...')
        handler.run(move=True)
    else:
        print('Nothing to install!')
        sys.exit(0)


def check_or_fail(handler):
    print(ui.header('Database Check', level=3))
    upofdate, missing = handler.check_uptodate()
    if outofdate:
        print(ui.paragraph('Must install databases to continue. To do so,'
                           ' run `dammit databases --install`. If you have'
                           ' already installed them, make sure you\'ve given'
                           ' the correct location to `--database-dir` or have'
                           ' exported the $DAMMIT_DB_DIR environment'
                           ' variable.'))
        sys.exit(2)


def register_builtin_tasks(handler, config, databases):

    settings = config['settings']

    register_pfam_tasks(handler, settings, databases)
    register_rfam_tasks(handler, settings, databases)
    register_orthodb_tasks(handler, settings, databases)
    register_busco_tasks(handler, settings, databases)
    register_uniref90_tasks(handler, settings, databases)

    return handler
                          

def register_pfam_tasks(handler, settings, databases):
    pfam_A = databases['Pfam-A']
    handler.register_task('download:Pfam-A',
                          get_download_and_gunzip_task(pfam_A['url'],
                                                       pfam_A['filename']),
                          files={'Pfam-A': pfam_A['filename']})
    handler.register_task('hmmpress:Pfam-A',
                          get_hmmpress_task(pfam_A['filename'], settings['hmmer']))
    return handler


def register_rfam_tasks(handler, settings, databases):
    rfam = databases['Rfam']
    handler.register_task('download:Rfam',
                           get_download_and_gunzip_task(rfam['url'], 
                                                        rfam['filename']),
                           files={'Rfam': rfam['filename']})
    handler.register_task('cmpress:Rfam',
                          get_cmpress_task(rfam['filename'], settings['infernal']))
    return handler


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
    return handler


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
    return handler


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
    return handler

