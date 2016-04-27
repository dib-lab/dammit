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
DOIT_VERBOSITY = 0


# Configuration stuff!
rel_path = os.path.dirname(__file__)
with open(os.path.join(rel_path, '.databases.json'), 'r') as fp:
    DATABASES = json.load(fp)
with open(os.path.join(rel_path, '.config.json'), 'r') as fp:
    CONFIG = json.load(fp)


def get_dammit_dir():
    return os.path.join(os.environ['HOME'],
                        CONFIG['settings']['dammit_dir'])


def uniqify_tasks(tasks):
    seen = set()
    uniqified = []
    for task in tasks:
        if task.name not in seen:
            uniqified.append(task)
            seen.add(task.name)
    return uniqified


def run_tasks(tasks, args, config={'verbosity': 0}):
    
    if type(tasks) is not list:
        raise TypeError('tasks must be a list')

    tasks = uniqify_tasks(tasks)
   
    class Loader(TaskLoader):
        @staticmethod
        def load_tasks(cmd, opt_values, pos_args):
            return tasks, config
   
    return DoitMain(Loader()).run(args)




