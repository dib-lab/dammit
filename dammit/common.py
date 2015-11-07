#!/usr/bin/env python
from __future__ import print_function

import json
import logging
import os
import sys
import textwrap
import time

from doit.cmd_base import TaskLoader
from doit.doit_cmd import DoitMain

CUR_TIME = time.strftime('%Y-%m-%d-%H%M')

DOIT_BACKEND = 'sqlite3'
DOIT_VERBOSITY = 2


# Configuration stuff!
rel_path = os.path.dirname(__file__)
with open(os.path.join(rel_path, '.databases.json'), 'r') as fp:
    DATABASES = json.load(fp)
with open(os.path.join(rel_path, '.config.json'), 'r') as fp:
    CONFIG = json.load(fp)


def get_dammit_dir():
    return os.path.join(os.environ['HOME'],
                        CONFIG['settings']['dammit_dir'])


def run_tasks(tasks, args, config={'verbosity': 0}):
    
    if type(tasks) is not list:
        raise TypeError('tasks must be a list')
   
    class Loader(TaskLoader):
        @staticmethod
        def load_tasks(cmd, opt_values, pos_args):
            return tasks, config
   
    DoitMain(Loader()).run(args)


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


def add_header(msg, level):
    s = ''
    if level == 0:
        symbol = '='
        return '{0}\n{1}\n{2}\n'.format(symbol * 40,
                                      msg,
                                      symbol * 40)
    elif level == 1:
        symbol = '~'
        return '{0}\n{1}\n{2}'.format(symbol * 30,
                                      msg,
                                      symbol * 30)
    else:
        symbol = '---'
        return '\n{0} {1}\n'.format(symbol,
                                    msg)


def print_header(msg, level):
    '''Standardize output headers for submodules.

    This doesn't need to be logged, but it's nice for
    the user.
    '''
    print(add_header(msg, level), file=sys.stderr)



