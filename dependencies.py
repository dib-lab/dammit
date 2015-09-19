#!/usr/bin/env python
from __future__ import print_function

from platform import system
from os.path import exists
import sys

from tasks import *

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

def get_dammit_dir(config):
    return os.path.join(os.eviron['HOME'],
                        config['settings']['dammit_dir'])


def check_dependencies(paths):
    '''Check that the given paths exist.

    
    '''
    for key, filename in paths:
        if not exists(filename):
            raise RuntimeError(
                '*** Failed to confirm path for {key}: {filename}. '\
                'Has the config file been modified?')


def get_dependencies_tasks(config):
    '''Generate tasks for installing the dependencies.

    These tasks download the dependencies and unpack them. Unlike the database
    tasks, these do not have the option of being put anywhere but in
    $HOME/.dammit/dependencies. The current binary dependencies included with
    dammit are:

        * BUSCO v1.1b
        * HMMER 3.1b2
        * Infernal 1.1.1
        * BLAST+ 2.2.31
        * TransDecoder 2.0.1

    This function also detects the operating system and selects the correct
    paths accordingly.

    Args:
        config (dict): The config dictionary.

    Returns:
    dict: A dictionary of the final binary paths.
    list: The doit tasks.

    '''

    dep_dir = os.path.join(get_dammit_dir(), config['settings']['dep_dir'])

    # This fails hard on Windows. For shame.
    # It also might fail hard on some linux distros, so probably should
    # actually looked at more closely.
    cur_platform = 'macosx' if system() == 'Darwin' else 'linux'
    
    tasks = []
    paths = {}

    tasks.append(
        get_download_and_untar_task(config['settings']['dep_url'],
                                    dep_dir,
                                    label='dependencies')
    )

    try:
        paths['HMMER'] = config['settings']['hmmer']['path'][cur_platform]
        paths['INFERNAL'] = config['settings']['infernal']['path'][cur_platform]
        paths['BUSCO'] = config['settings']['busco']['path']
        paths['TRANSDECODER'] = config['settings']['transdecoder']['path']
        paths['BLAST'] = 
    except KeyError as e:
        raise RuntimeError(
            'Something went wrong accessing the config dictionary --'\
            ' has it been modified?')

    return paths, tasks

