#!/usr/bin/enb python
from __future__ import print_function

import logging
import os
from os import path
import sys

from doit.dependency import Dependency, SqliteDB

from . import ui
from .handler import TaskHandler
from .log import LogReporter
from .hmmer import get_hmmpress_task
from .infernal import get_cmpress_task
from .last import get_lastdb_task
from .tasks import get_download_and_gunzip_task, \
                   get_download_and_untar_task

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

    handler = TaskHandler(directory, logger, config=config,
                          backend=config['doit_backend'],
                          verbosity=args.verbosity)

    return register_builtin_tasks(handler, config, databases,
                                  with_uniref=args.full)


def install(handler):
    uptodate, missing = handler.print_uptodate()
    if not uptodate:
        print('Installing...')
        handler.run(move=True)
    else:
        print('Nothing to install!')
        sys.exit(0)


def check_or_fail(handler):
    print(ui.header('Database Check', level=3))
    uptodate, missing = handler.print_uptodate()
    if not uptodate:
        print(ui.paragraph('Must install databases to continue. To do so,'
                           ' run `dammit databases --install`. If you have'
                           ' already installed them, make sure you\'ve given'
                           ' the correct location to `--database-dir` or have'
                           ' exported the $DAMMIT_DB_DIR environment'
                           ' variable.'))
        sys.exit(2)


def register_builtin_tasks(handler, config, databases, with_uniref=False):

    register_pfam_tasks(handler, config['hmmer']['hmmpress'], databases)
    register_rfam_tasks(handler, config['infernal']['cmpress'], databases)
    register_orthodb_tasks(handler, config['last']['lastdb'], databases)
    register_busco_tasks(handler, config['busco'], databases)
    if with_uniref:
        register_uniref90_tasks(handler, config['last']['lastdb'], databases)

    return handler
                          

def register_pfam_tasks(handler, params, databases):
    pfam_A = databases['Pfam-A']
    task = get_download_and_gunzip_task(pfam_A['url'],
                                        pfam_A['filename'])
    handler.register_task('download:Pfam-A',
                          task,
                          files={'Pfam-A': pfam_A['filename']})
    handler.register_task('hmmpress:Pfam-A',
                          get_hmmpress_task(pfam_A['filename'], 
                                            params=params,
                                            task_dep=[task.name]))
    return handler


def register_rfam_tasks(handler, params, databases):
    rfam = databases['Rfam']
    task = get_download_and_gunzip_task(rfam['url'], 
                                        rfam['filename'])
    handler.register_task('download:Rfam',
                           task,
                           files={'Rfam': rfam['filename']})
    handler.register_task('cmpress:Rfam',
                          get_cmpress_task(rfam['filename'],
                                           task_dep=[task.name],
                                           params=params))
    return handler


def register_orthodb_tasks(handler, params, databases):
    orthodb = databases['OrthoDB']
    task = get_download_and_gunzip_task(orthodb['url'], 
                                        orthodb['filename'])
    handler.register_task('download:OrthoDB',
                          task,
                          files={'OrthoDB': orthodb['filename']})
    handler.register_task('lastdb:OrthoDB',
                          get_lastdb_task(orthodb['filename'], 
                                          orthodb['filename'], 
                                          prot=True,
                                          params=params,
                                          task_dep=[task.name]))
    return handler


def register_busco_tasks(handler, config, databases):
    busco = databases['BUSCO']
    busco_dir = config['db_dir']
    for group_name in busco:
        group = busco[group_name]
        files = {'BUSCO-{0}'.format(group_name): path.join(busco_dir, group_name)}
        handler.register_task('download:BUSCO-{0}'.format(group_name),
                              get_download_and_untar_task(group['url'],
                                                         busco_dir,
                                                         label=group_name),
                              files=files)
    return handler


def register_uniref90_tasks(handler, params, databases):
    uniref90 = databases['uniref90']
    task = get_download_and_gunzip_task(uniref90['url'],
                                        uniref90['filename'])
    handler.register_task('download:uniref90',
                          task,
                          files={'uniref90': uniref90['filename']})
    handler.register_task('lastdb:uniref90',
                          get_lastdb_task(uniref['filename'],
                                          uniref['filename'],
                                          prot=True,
                                          params=params,
                                          task_dep=[task.name]))
    return handler

