#!/usr/bin/enb python
from __future__ import print_function

import os

import common
from tasks import get_download_and_gunzip_task, \
                  get_hmmpress_task, \
                  get_cmpress_task, \
                  get_blast_format_task, \
                  get_download_and_untar_task, \
                  get_lastdb_task


def get_database_tasks(db_dir, prog_paths, busco_db, full):
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

    try:
        os.mkdir(db_dir)
    except OSError:
        pass

    tasks = []
    databases = {}

    # Get Pfam-A and prepare it for use with hmmer
    PFAM = os.path.join(db_dir, common.DATABASES['pfam']['filename'])
    tasks.append(
        get_download_and_gunzip_task(common.DATABASES['pfam']['url'], PFAM)
    )
    hmmer_dir = prog_paths.get('hmmer', '')
    tasks.append(
        get_hmmpress_task(PFAM, common.CONFIG['settings']['hmmer'],
                          hmmer_dir=hmmer_dir)
    )
    databases['PFAM'] = os.path.abspath(PFAM)

    # Get Rfam and prepare it for use with Infernal
    RFAM = os.path.join(db_dir, common.DATABASES['rfam']['filename'])
    infernal_dir = prog_paths.get('infernal', '')
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
    '''
    blast_dir = prog_paths.get('blast', '')
    tasks.append(
        get_blast_format_task(ORTHODB, ORTHODB + '.db', 
                              common.DATABASES['orthodb']['db_type'],
                              blast_dir=blast_dir)
    )'''
    last_dir = prog_paths.get('last', '')
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

def run_database_tasks(db_dir, tasks, args=['run']):
    
    doit_config = {
                    'backend': common.DOIT_BACKEND,
                    'verbosity': common.DOIT_VERBOSITY,
                    'dep_file': os.path.join(db_dir, 'databases.doit.db')
                  }

    common.run_tasks(tasks, args, config=doit_config)


