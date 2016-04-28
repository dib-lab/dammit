#!/usr/bin/env python

import os

from doit.cmd_base import TaskLoader
from doit.doit_cmd import DoitMain
from doit.dependency import Dependency, SqliteDB

from .utils import Move
from . import ui

class TaskHandler(TaskLoader):

    def __init__(self, directory, logger, config=None, files=None, **doit_config_kwds):

        super(TaskHandler, self).__init__()

        if files is None:
            self.files = {}
        elif type(files) is not dict:
            raise TypeError('files must be of type dict')
        else:
            self.files = files

        self.tasks = {}

        self.directory = os.path.abspath(directory)

        try:
            os.mkdir(self.directory)
        except OSError:
            logger.debug('database dir already exists')

        dep_file = os.path.join(self.directory, 'doit.db')
        self.doit_config = dict(dep_file=dep_file, **doit_config_kwds)
        self.doit_dep_mgr = Dependency(SqliteDB, dep_file)
        self.logger = logger

    def register_task(self, name, task, files=None):
        if files is None:
            files = {}
        if type(files) is not dict:
            raise TypeError('files must be of type dict')
        
        self.tasks[name] = task
        self.files.update(files)
        self.logger.debug('registered task {0}: {1}\n'
                          '  with files {2}'.format(name, task, files))

    def clear_tasks(self):
        self.logger.debug('Clearing {0} tasks'.format(len(self.tasks)))
        self.tasks = {}

    def get_status(self, task):
        if type(task) is str:
            try:
                task = self.tasks[task]
            except KeyError:
                self.logger.error('Task not found:{0}'.format(task))
                raise
        status = self.doit_dep_mgr.get_status(task, self.tasks.values()).status
        self.logger.debug('Task {0} had status {1}'.format(task, status))
        return status

    def print_uptodate(self):
        uptodate, outofdate = self.check_uptodate()
        if uptodate:
            print(ui.paragrap('All tasks up-to-date!'))
        else:
            print('\nOut-of-date tasks:')
            print(ui.listing(outofdate))
        return uptodate, outofdate

    def check_uptodate(self):
        outofdate_tasks = {}
        outofdate = False
        for task_name, task in self.tasks.items():
            status = self.get_status(task)
            if status != 'up-to-date':
                outofdate_tasks[task_name] = status
                outofdate = True
        return not outofdate, outofdate_tasks
    
    def load_tasks(self, cmd, opt_values, pos_args):
        self.logger.debug('loading {0} tasks'.format(len(self.tasks)))
        return self.tasks.values()

    def run(self, doit_args=None, move=False):
        if doit_args is None:
            doit_args = ['run']
        runner = DoitMain(self())
        if move:
            with Move(self.directory):
                return runner.run(doit_args)
        else:
            return runner.run(doit_args)
