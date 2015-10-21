#!/usr/bin/enb python
from __future__ import print_function

import logging
import os

from doit.dependency import Dependency, SqliteDB

import common
from tasks import get_download_and_gunzip_task, \
                  get_hmmpress_task, \
                  get_cmpress_task, \
                  get_download_and_untar_task, \
                  get_lastdb_task

logger = logging.getLogger(__name__)

def get_dir(args):
    # By default, we store databases in the home directory
    db_dir = args.database_dir
    if db_dir is None:
        db_dir = os.path.join(common.get_dammit_dir(), 
                              common.CONFIG['settings']['db_dir'])
    return os.path.abspath(db_dir)

def get_tasks(db_dir, args):
    '''Generate tasks for installing the bundled databases. 
    
    These tasks download the databases, unpack them, and format them for use.
    Current bundled databases are:

        * Pfam-A (protein domans)
        * Rfam (RNA models)
        * OrthoDB8 (conserved ortholog groups)
        * uniref90 (protiens, if --full selected)
    
    User-supplied databases are downloaded separately.

    Args:
        db_dir (str): Directory where the databases will be stored.
        busco_db (str): The BUSCO group to use.
        full (bool): Whether to do a full run and get UNIREF90 as well.

    Returns:
        dict: A dictionary of the final database paths.
        list: A list of the doit tasks.

    '''

    tasks = []
    databases = {}

    # Get Pfam-A and prepare it for use with hmmer
    PFAM = os.path.join(db_dir, common.DATABASES['pfam']['filename'])
    tasks.append(
        get_download_and_gunzip_task(common.DATABASES['pfam']['url'], PFAM)
    )
    tasks.append(
        get_hmmpress_task(PFAM, common.CONFIG['settings']['hmmer'])
    )
    databases['PFAM'] = os.path.abspath(PFAM)

    # Get Rfam and prepare it for use with Infernal
    RFAM = os.path.join(db_dir, common.DATABASES['rfam']['filename'])
    tasks.append(
        get_download_and_gunzip_task(common.DATABASES['rfam']['url'], RFAM)
    )
    tasks.append(
        get_cmpress_task(RFAM, common.CONFIG['settings']['infernal'])
    )
    databases['RFAM'] = os.path.abspath(RFAM)

    # Get OrthoDB and prepare it for BLAST use
    ORTHODB = os.path.join(db_dir, common.DATABASES['orthodb']['filename'])
    tasks.append(
        get_download_and_gunzip_task(common.DATABASES['orthodb']['url'], ORTHODB)
    )

    lastdb_cfg = common.CONFIG['settings']['last']['lastdb']
    tasks.append(
        get_lastdb_task(ORTHODB, ORTHODB + '.db', lastdb_cfg, prot=True)
    )
    ORTHODB += '.db'
    databases['ORTHODB'] = os.path.abspath(ORTHODB)

    ORTHODB_GENES = os.path.join(db_dir,
                                 common.DATABASES['orthodb_genes']['filename'])
    tasks.append(
        get_download_and_gunzip_task(common.DATABASES['orthodb_genes']['url'],
                                     ORTHODB_GENES)
    )
    databases['ORTHODB_GENES'] = os.path.abspath(ORTHODB_GENES)

    # A little confusing. First, we get the top-level BUSCO path:
    BUSCO = os.path.join(db_dir, 'buscodb')
    tasks.append(
        # That top-level path is given to the download task:
        get_download_and_untar_task(common.DATABASES['busco'][args.busco_group]['url'], 
                                    BUSCO,
                                    label=args.busco_group)
    )
    # The untarred arhive has a folder named after the group:
    databases['BUSCO'] = os.path.abspath(os.path.join(BUSCO, args.busco_group))

    # Get uniref90 if the user specifies
    if args.full:
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

def get_doit_config(db_dir):

    doit_config = {
                    'reporter': common.LogReporter(logger),
                    'backend': common.DOIT_BACKEND,
                    'verbosity': common.DOIT_VERBOSITY,
                    'dep_file': os.path.join(db_dir, 'databases.doit.db')
                  }
    logger.debug('doit_config:{0}'.format(doit_config))
    return doit_config

def run_tasks(db_dir, tasks, args=['run']):
    
    doit_config = get_doit_config(db_dir)
    common.run_tasks(tasks, args, config=doit_config)

def do_check(db_dir, args):

    common.print_header('Checking for database prep (dir: {0})'.format(db_dir),
                        level=2)

    databases, tasks = get_tasks(db_dir, args)

    doit_config = get_doit_config(db_dir)
    dep_manager = Dependency(SqliteDB, doit_config['dep_file'])
    missing = False
    for task in tasks:
        status = dep_manager.get_status(task, tasks)
        logger.debug('{0}:{1}'.format(task.name, status.status))
        if status.status != 'up-to-date':
            missing = True
            logger.warning('[ ] {0}'.format(task.name))
        else:
            logger.info('[x] {0}'.format(task.name))

    common.print_header('Database results', level=2)

    if missing:
        logger.warning('Database prep incomplete')
        common.print_header('to prepare databases, run: dammit databases'\
                            ' --install', level=2)
    else:
        logger.info('All databases prepared!')

    return databases, tasks, missing
