#!/usr/bin/env python
from __future__ import print_function

import os
from platform import system
import sys

import common
from tasks import get_download_and_untar_task

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

def get_dammit_dir():
    return os.path.join(os.environ['HOME'],
                        common.CONFIG['settings']['dammit_dir'])

def get_dependency_dir():
    return os.path.join(get_dammit_dir(), \
                        common.CONFIG['settings']['dep_dir'])

def get_dependency_tasks():
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

    dep_dir = get_dependency_dir()

    # This fails hard on Windows. For shame.
    # It also might fail hard on some linux distros, so probably should
    # actually looked at more closely.
    cur_platform = 'macosx' if system() == 'Darwin' else 'linux'
    
    tasks = []
    paths = {}

    # Striping the .tar.gz
    strip_ext = lambda s: s[:-7]

    # Check for hmmer
    hmmscan = which('hmmscan')
    hmmpress = which('hmmpress')
    if hmmscan is None or hmmpress is None:
        url = common.CONFIG['settings']['hmmer']['url'][cur_platform]
        tasks.append(
            get_download_and_untar_task(url,
                                        dep_dir,
                                        label='hmmer')
        )
        hmmer_dir = strip_ext(os.path.basename(url))
        paths['hmmer'] = os.path.join(dep_dir, hmmer_dir)

    # Check for infernal
    cmscan = which('cmscan')
    cmpress = which('cmpress')
    if cmscan is None or cmpress is None:
        url = common.CONFIG['settings']['infernal']['url'][cur_platform]
        tasks.append(
            get_download_and_untar_task(url,
                                        dep_dir,
                                        label='infernal')
        )
        infernal_dir = strip_ext(os.path.basename(url))
        paths['infernal'] = os.path.join(dep_dir, infernal_dir)

    blastn = which('blastn')
    blastx = which('blastx')
    tblastn = which('tblastn')
    makeblastdb = which('makeblastdb')
    if blastn is None or blastx is None \
        or tblastn is None or makeblastdb is None:

        url = common.CONFIG['settings']['blast']['url'][cur_platform]
        tasks.append(
            get_download_and_untar_task(url,
                                        dep_dir,
                                        label='blast')
        )
        blast_dir = strip_ext(os.path.basename(url))
        paths['blast'] = os.path.join(dep_dir, blast_dir)

    # Check for BUSCO
    busco = which('BUSCO_v1.1b1.py')
    if busco is None:
        url = common.CONFIG['settings']['busco']['url']
        tasks.append(
            get_download_and_untar_task(url,
                                        dep_dir,
                                        label='busco')
        )
        busco_dir = strip_ext(os.path.basename(url))
        paths['busco'] = os.path.join(dep_dir, busco_dir)

    longorfs = which('TransDecoder.LongOrfs')
    predict = which('TransDecoder.Predict')
    if longorfs is None or predict is None:
        url = common.CONFIG['settings']['busco']['url']
        tasks.append(
            get_download_and_untar_task(url,
                                        dep_dir,
                                        label='transdecoder')
        )
        transdecoder_dir = strip_ext(os.path.basename(url))
        paths['transdecoder'] = os.path.join(dep_dir, transdecoder_dir)

    return paths, tasks


def run_dependency_tasks(tasks, args=['run']):
    '''
    This set of tasks keeps its own doit db in the db folder to share
    them between all runs.
    '''
    
    dep_dir = get_dependency_dir()
    doit_config = {
                    'backend': common.DOIT_BACKEND,
                    'verbosity': common.DOIT_VERBOSITY,
                    'dep_file': os.path.join(dep_dir, 'dependencies.doit.db')
                  }

    common.run_tasks(tasks, args, config=doit_config)
