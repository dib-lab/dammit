#!/usr/bin/env python
from __future__ import print_function

import json
import logging
import os
import sys
import time

from doit.cmd_base import TaskLoader
from doit.doit_cmd import DoitMain

CUR_TIME = time.strftime('%Y-%m-%d-%H%M')

DOIT_VERBOSITY = 2

# Configuration stuff!
with open('.databases.json', 'r') as fp:
    DATABASES = json.load(fp)
with open('.config.json', 'r') as fp:
    CONFIG = json.load(fp)

def get_dammit_dir():
    return os.path.join(os.environ['HOME'],
                        CONFIG['settings']['dammit_dir'])

def run_tasks(tasks, args, config={'verbosity': 2}):
    
    if type(tasks) is not list:
        raise TypeError('tasks must be a list')
   
    class Loader(TaskLoader):
        @staticmethod
        def load_tasks(cmd, opt_values, pos_args):
            return tasks, config
   
    DoitMain(Loader()).run(args)

def print_header(msg, level):
    '''Standardize output headers for submodules.

    This doesn't need to be logged, but it's nice for
    the user.
    '''
    if level == 0:
        symbol = '='
        N = 40
    elif level == 1:
        symbol = '~'
        N = 30
    else:
        symbol = '-'
        N = 20

    print(symbol * N, file=sys.stderr)
    print(msg, file=sys.stderr)
    print(symbol * N, '\n', file=sys.stderr)

log_dir = os.path.join(get_dammit_dir(), 'log')
try:
    os.makedirs(log_dir)
except OSError:
    pass
log_file = os.path.join(log_dir, 'dammit-all.log')

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(name)s:%(funcName)s:%(lineno)d '\
                           '[%(levelname)s] \n%(message)s\n-----',
                    datefmt='%m-%d %H:%M:%S',
                    filename=log_file,
                    filemode='a')

console = logging.StreamHandler(sys.stderr)
console.setLevel(logging.INFO)
formatter = logging.Formatter('\t%(message)-12s [%(name)s:%(levelname)s]')
console.setFormatter(formatter)
logging.getLogger('').addHandler(console)

logging.getLogger('').debug('*** dammit! begin ***')
