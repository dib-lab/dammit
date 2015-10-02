#!/usr/bin/env python
from __future__ import print_function

import logging
import os
from platform import system
import sys

from doit.dependency import Dependency, SqliteDB

import common
from tasks import get_download_and_untar_task

logger = logging.getLogger(__name__)

def which(program):
    '''Checks whether the given program (or program path) is valid and
    executable.

    NOTE: Sometimes copypasta is okay! This function came from stackoverflow:

        http://stackoverflow.com/a/377028/5109965

    Args:
        program (str): Either a program name or full path to a program.

    Returns:
        Return the path to the executable or None if not found
    '''
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None

def check_system():
    
    deps = {}

    hmmscan = which('hmmscan')
    hmmpress = which('hmmpress')
    if hmmscan is None or hmmpress is None:
        deps['HMMER'] = False
    else:
        deps['HMMER'] = True

    cmscan = which('cmscan')
    cmpress = which('cmpress')
    if cmscan is None or cmpress is None:
        deps['Infernal'] = False
    else:
        deps['Infernal'] = True

    blastn = which('blastn')
    blastx = which('blastx')
    tblastn = which('tblastn')
    makeblastdb = which('makeblastdb')
    if blastn is None or blastx is None \
        or tblastn is None or makeblastdb is None:

        deps['BLAST+'] = False
    else:
        deps['BLAST+'] = True

    busco = which('BUSCO_v1.1b1.py')
    if busco is None:
        deps['BUSCO'] = False
    else:
        deps['BUSCO'] = True


    longorfs = which('TransDecoder.LongOrfs')
    predict = which('TransDecoder.Predict')
    if longorfs is None or predict is None:
        deps['TransDecoder'] = False
    else:
        deps['TransDecoder'] = True

    lastdb = which('lastdb')
    lastal = which('lastal')
    if lastdb is None or lastal is None:
        deps['LAST'] = False
    else:
        deps['LAST'] = True

    return deps

def get_dir():
    return os.path.join(common.get_dammit_dir(), \
                        common.CONFIG['settings']['dep_dir'])

def get_tasks():
    '''Check for dependencies and generate tasks for them if missing.

    These tasks check for each dependency on the system PATH, and generate a
    task to download and unpack them if they're missing. Unlike the database
    tasks, these do not have the option of being put anywhere but in
    $HOME/.dammit/dependencies. The current binary dependencencies which dammit
    will download when missing are:

        * BUSCO v1.1b
        * HMMER 3.1b2
        * Infernal 1.1.1
        * BLAST+ 2.2.31
        * TransDecoder 2.0.1

    Returns:
    dict: A dictionary of the final binary paths.
    list: The doit tasks.

    '''

    dep_dir = get_dir()

    # This fails hard on Windows. For shame.
    # It also might fail hard on some linux distros, so probably should
    # actually looked at more closely.
    cur_platform = 'macosx' if system() == 'Darwin' else 'linux'
    
    system_deps = check_system()
    dammit_deps = {}
    tasks = {}

    # Striping the .tar.gz
    strip_ext = lambda s: s[:-7]

    # Check for hmmer
    if not system_deps['HMMER']:
        url = common.CONFIG['settings']['hmmer']['url'][cur_platform]
        tasks['HMMER'] = get_download_and_untar_task(url,
                                                     dep_dir,
                                                     label='hmmer')
        hmmer_dir = strip_ext(os.path.basename(url))
        dammit_deps['HMMER'] = os.path.join(dep_dir, hmmer_dir)

    # Check for infernal
    if not system_deps['Infernal']:
        url = common.CONFIG['settings']['infernal']['url'][cur_platform]
        tasks['Infernal'] = get_download_and_untar_task(url,
                                                        dep_dir,
                                                        label='infernal')
        infernal_dir = strip_ext(os.path.basename(url))
        dammit_deps['Infernal'] = os.path.join(dep_dir, infernal_dir)

    if not system_deps['BLAST+']:
        url = common.CONFIG['settings']['blast']['url'][cur_platform]
        tasks['BLAST+'] = get_download_and_untar_task(url,
                                                      dep_dir,
                                                      label='blast')
        blast_dir = strip_ext(os.path.basename(url))
        dammit_deps['BLAST+'] = os.path.join(dep_dir, blast_dir)

    # Check for BUSCO
    if not system_deps['BUSCO']:
        url = common.CONFIG['settings']['busco']['url']
        tasks['BUSCO'] = get_download_and_untar_task(url,
                                                     dep_dir,
                                                     label='busco')
        busco_dir = strip_ext(os.path.basename(url))
        dammit_deps['BUSCO'] = os.path.join(dep_dir, busco_dir)

    if not system_deps['TransDecoder']:
        url = common.CONFIG['settings']['transdecoder']['url']
        tasks['TransDecoder'] = get_download_and_untar_task(url,
                                                            dep_dir,
                                                            label='transdecoder')
        transdecoder_dir = strip_ext(os.path.basename(url))
        dammit_deps['TransDecoder'] = os.path.join(dep_dir, transdecoder_dir)

    if not system_deps['LAST']:
        url = common.CONFIG['settings']['last']['url'][cur_platform]
        tasks['LAST'] = get_download_and_untar_task(url,
                                                    dep_dir,
                                                    label='last')
        last_dir = strip_ext(os.path.basename(url))
        dammit_deps['LAST'] = os.path.join(dep_dir, last_dir)

    return dammit_deps, system_deps, tasks

def get_doit_config():
    
    dep_dir = get_dir()
    doit_config = {
                    'backend': common.DOIT_BACKEND,
                    'verbosity': common.DOIT_VERBOSITY,
                    'dep_file': os.path.join(dep_dir, 'dependencies.doit.db')
                  }

    return doit_config


def run_tasks(tasks, args=['run']):
    '''
    This set of tasks keeps its own doit db in the db folder to share
    them between all runs.
    '''

    doit_config = get_doit_config()
    common.run_tasks(tasks, args, config=doit_config)

def check():

    dammit_deps, system_deps, dep_tasks = get_tasks()

    dammit_dep_status = {}

    doit_config = get_doit_config()
    dep_manager = Dependency(SqliteDB, doit_config['dep_file'])
    for key in system_deps:
        try:
            status = dep_manager.get_status(dep_tasks[key], dep_tasks)
        except KeyError:
            dammit_dep_status[key] = False
        else:
            if status != 'up-to-date':
                dammit_dep_status[key] = False
            else:
                dammit_dep_status[key] = True

    missing = []
    for ((key, system_status), (_, dammit_status)) in zip(system_deps.items(),
                                                          dammit_dep_status.items()):
        if system_status == dammit_status and not system_status:
            missing.append(key)
        elif dammit_status:
            logger.info('{0} found [previously installed by dammit]'.format(key))
        else:
            logger.info('{0} found [installed on system PATH]'.format(key))

    if missing:
        logger.warning('* {0} missing'.format(', '.join(missing)))
        common.print_header('to get dependencies, run: dammit dependencies --install',
                            level=2)
    else:
        logger.info('* all dependencies satisfied!')

    return missing

