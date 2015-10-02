#!/usr/bin/enb python
from __future__ import print_function

import logging
import os

from doit.dependency import Dependency, SqliteDB

import common
from tasks import get_download_and_gunzip_task, \
                  get_hmmpress_task, \
                  get_cmpress_task, \
                  get_blast_format_task, \
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

def get_tasks(db_dir, prog_paths, busco_db, full):
    '''Generate tasks for installing the bundled databases. 
    
    These tasks download the databases, unpack them, and format them for use.
    Current bundled databases are:

        * Pfam-A (protein domans)
        * Rfam (RNA models)
        * OrthoDB8 (conserved ortholog groups)
        * uniref90 (protiens, if --full selected)
    
    User-supplied databases are downloaded separately.

    Args:
        db_dir (str): Directory where the databases wille be stored.
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
    hmmer_dir = prog_paths.get('HMMER', '')
    tasks.append(
        get_hmmpress_task(PFAM, common.CONFIG['settings']['hmmer'],
                          hmmer_dir=hmmer_dir)
    )
    databases['PFAM'] = os.path.abspath(PFAM)

    # Get Rfam and prepare it for use with Infernal
    RFAM = os.path.join(db_dir, common.DATABASES['rfam']['filename'])
    infernal_dir = prog_paths.get('Infernal', '')
    tasks.append(
        get_download_and_gunzip_task(common.DATABASES['rfam']['url'], RFAM)
    )
    tasks.append(
        get_cmpress_task(RFAM, common.CONFIG['settings']['infernal'],
                         infernal_dir=infernal_dir)
    )
    databases['RFAM'] = os.path.abspath(RFAM)

    # Get OrthoDB and prepare it for BLAST use
    ORTHODB = os.path.join(db_dir, common.DATABASES['orthodb']['filename'])
    tasks.append(
        get_download_and_gunzip_task(common.DATABASES['orthodb']['url'], ORTHODB)
    )

    blast_dir = prog_paths.get('BLAST+', '')

    last_dir = prog_paths.get('LAST', '')
    lastdb_cfg = common.CONFIG['settings']['last']['lastdb']
    tasks.append(
        get_lastdb_task(ORTHODB, ORTHODB + '.db', lastdb_cfg,
                        prot=True, last_dir=last_dir)
    )
    ORTHODB += '.db'
    databases['ORTHODB'] = os.path.abspath(ORTHODB)

    # A little confusing. First, we get the top-level BUSCO path:
    BUSCO = os.path.join(db_dir, 'buscodb')
    tasks.append(
        # That top-level path is given to the download task:
        get_download_and_untar_task(common.DATABASES['busco'][busco_db]['url'], 
                                    BUSCO,
                                    label=busco_db)
    )
    # The untarred arhive has a folder named after the group:
    databases['BUSCO'] = os.path.abspath(os.path.join(BUSCO, busco_db))

    # Get uniref90 if the user specifies
    if full:
        UNIREF = os.path.join(common.DATABASES['uniref90']['filename'])
        tasks.append(
            get_download_and_gunzip_task(common.DATABASES['uniref90']['url'], UNIREF)
        )
        tasks.append(
            get_blast_format_task(UNIREF, UNIREF + '.db',
                                  common.DATABASES['uniref90']['db_type'],
                                  blast_dir=blast_dir)
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
    return doit_config

def run_tasks(db_dir, tasks, args=['run']):
    
    doit_config = get_doit_config(db_dir)
    common.run_tasks(tasks, args, config=doit_config)

def check(tasks, db_dir):

    doit_config = get_doit_config(db_dir)
    dep_manager = Dependency(SqliteDB, doit_config['dep_file'])
    missing = False
    for task in tasks:
        status = dep_manager.get_status(task, tasks)
        if status != 'up-to-date':
            missing = True
            logger.warning('{0} not ready'.format(task.name))
        else:
            logger.info('{0} ready!'.format(task.name))

    if missing:
        logger.warning('* database prep incomplete')
        common.print_header('to prepare databases, run: dammit databases'\
                            ' --install', level=2)
    else:
        logger.info('* all databases prepared!')

    return missing
