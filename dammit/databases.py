#!/usr/bin/enb python
from __future__ import print_function

import logging
import os
from os import path
import sys

from doit.dependency import Dependency, SqliteDB
from shmlast.last import lastdb_task as get_lastdb_task

from . import ui
from .handler import TaskHandler
from .tasks.hmmer import get_hmmpress_task
from .tasks.infernal import get_cmpress_task
from .tasks.shell import (get_download_and_gunzip_task,
                          get_download_and_untar_task)

def get_handler(args, config, databases):

    logger = logging.getLogger('DatabaseHandler')
    logger.debug('get_handler')

    handler = TaskHandler(args.database_dir, logger, config=config,
                          db='databases',
                          backend=config['doit_backend'],
                          verbosity=args.verbosity)

    return register_builtin_tasks(handler, config, databases,
                                  with_uniref=args.full)


def default_database_dir(logger):
    try:
        directory = os.environ['DAMMIT_DB_DIR']
        logger.debug('found DAMMIT_DB_DIR env variable')
    except KeyError:
        logger.debug('no DAMMIT_DB_DIR or --database-dir, using'\
                          ' default')
        directory = path.join(os.environ['HOME'], '.dammit', 'databases')
    return directory


def install(handler):
    print(ui.header('Database Install', level=3))
    uptodate, missing = handler.print_uptodate()
    if not uptodate:
        print('Installing...')
        handler.run(move=True)
    else:
        print('Nothing to install!')
        sys.exit(0)


def check_or_fail(handler):
    print(ui.header('Database Check', level=3))
    print('Doit Database: {0}'.format(handler.dep_file))
    print('Database Directory: {0}'.format(handler.directory))
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
    filename = path.join(handler.directory, pfam_A['filename'])
    task = get_download_and_gunzip_task(pfam_A['url'],
                                        filename)
    handler.register_task('download:Pfam-A', task,
                          files={'Pfam-A': filename})
    handler.register_task('hmmpress:Pfam-A',
                          get_hmmpress_task(filename, 
                                            params=params,
                                            task_dep=[task.name]))
    return handler


def register_rfam_tasks(handler, params, databases):
    rfam = databases['Rfam']
    filename = path.join(handler.directory, rfam['filename'])
    task = get_download_and_gunzip_task(rfam['url'], 
                                        filename)
    handler.register_task('download:Rfam', task,
                           files={'Rfam': filename})
    handler.register_task('cmpress:Rfam',
                          get_cmpress_task(filename,
                                           task_dep=[task.name],
                                           params=params))
    return handler


def register_orthodb_tasks(handler, params, databases):
    orthodb = databases['OrthoDB']
    filename = path.join(handler.directory, orthodb['filename'])
    task = get_download_and_gunzip_task(orthodb['url'], 
                                        filename)
    handler.register_task('download:OrthoDB', task,
                          files={'OrthoDB': filename})
    handler.register_task('lastdb:OrthoDB',
                          get_lastdb_task(filename, 
                                          filename, 
                                          prot=True,
                                          params=params,
                                          task_dep=[task.name]))
    return handler


def register_busco_tasks(handler, config, databases):
    busco = databases['BUSCO']
    busco_dir = path.join(handler.directory, config['db_dir'])
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
    filename = path.join(handler.directory, uniref90['filename'])
    handler.register_task('download:uniref90',
                          task,
                          files={'uniref90': filename})
    handler.register_task('lastdb:uniref90',
                          get_lastdb_task(filename,
                                          filename,
                                          prot=True,
                                          params=params,
                                          task_dep=[task.name]))
    return handler

