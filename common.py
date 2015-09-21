#!/usr/bin/env python
from __future__ import print_function

import json

from doit.cmd_base import TaskLoader
from doit.doit_cmd import DoitMain

DOIT_BACKEND = 'sqlite3'
DOIT_VERBOSITY = 2

# Configuration stuff!
with open('.databases.json', 'r') as fp:
    DATABASES = json.load(fp)
with open('.config.json', 'r') as fp:
    CONFIG = json.load(fp)

def run_tasks(tasks, args, config={'verbosity': 2}):
    
    if type(tasks) is not list:
        raise TypeError('tasks must be a list')
   
    class Loader(TaskLoader):
        @staticmethod
        def load_tasks(cmd, opt_values, pos_args):
            return tasks, config
   
    DoitMain(Loader()).run(args)
