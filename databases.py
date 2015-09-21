#!/usr/bin/enb python
from __future__ import print_function

from tasks import *

def get_databases_tasks(resources, db_dir, busco_db, full):
    '''Generate tasks for installing the bundled databases. 
    
    These tasks download the databases, unpack them, and format them for use.
    Current bundled databases are:

        * Pfam-A (protein domans)
        * Rfam (RNA models)
        * OrthoDB8 (conserved ortholog groups)
        * uniref90 (protiens, if --full selected)
    
    User-supplied databases are downloaded separately.

    Args:
        resources (dict): The resources dictionary.
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
    PFAM = os.path.join(db_dir, resources['pfam']['filename'])
    tasks.append(
        get_download_and_gunzip_task(resources['pfam']['url'], PFAM)
    )
    tasks.append(
        get_hmmpress_task(PFAM)
    )
    databases['PFAM'] = os.path.abspath(PFAM)

    # Get Rfam and prepare it for use with Infernal
    RFAM = os.path.join(db_dir, resources['rfam']['filename'])
    tasks.append(
        get_download_and_gunzip_task(resources['rfam']['url'], RFAM)
    )
    tasks.append(
        get_cmpress_task(RFAM)
    )
    databases['RFAM'] = os.path.abspath(RFAM)

    # Get OrthoDB and prepare it for BLAST use
    ORTHODB = os.path.join(db_dir, resources['orthodb']['filename'])
    tasks.append(
        get_download_and_gunzip_task(resources['orthodb']['url'], ORTHODB)
    )
    tasks.append(
        get_blast_format_task(ORTHODB, ORTHODB + '.db', 
                              resources['orthodb']['db_type'])
    )
    ORTHODB += '.db'
    databases['ORTHODB'] = os.path.abspath(ORTHODB)

    # A little confusing. First, we get the top-level BUSCO path:
    BUSCO = os.path.join(db_dir, 'buscodb')
    tasks.append(
        # That top-level path is given to the download task:
        get_download_and_untar_task(resources['busco'][busco_db]['url'], 
                                    BUSCO,
                                    label=busco_db)
    )
    # The untarred arhive has a folder named after the group:
    databases['BUSCO'] = os.path.abspath(os.path.join(BUSCO, busco_db))

    # Get uniref90 if the user specifies
    if full:
        UNIREF = os.path.join(resources['uniref90']['filename'])
        tasks.append(
            get_download_and_gunzip_task(resources['uniref90']['url'], UNIREF)
        )
        tasks.append(
            get_blast_format_task(UNIREF, UNIREF + '.db',
                                  resources['uniref90']['db_type'])
        )
        UNIREF += '.db'
        databases['UNIREF'] = os.path.abspath(UNIREF)

    return databases, tasks

def run_install_databases(db_dir, tasks, args=['run']):
    
    doit_config = {
                    'backend': DB_BACKEND,
                    'verbosity': DOIT_VERBOSITY,
                    'dep_file': os.path.join(db_dir, 'databases.doit.db')
                  }

    run_tasks(tasks, args, config=doit_config)


