# Copyright (C) 2015-2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the BSD license.  See the LICENSE file for details.

from collections import OrderedDict
from os import mkdir, path

from doit.cmd_base import TaskLoader
from doit.doit_cmd import DoitMain
from doit.dependency import Dependency, SqliteDB

from dammit.profile import StartProfiler
from dammit.utils import Move
from dammit import ui


class TaskHandler(TaskLoader):

    def __init__(self, directory, logger, files=None, 
                 profile=False, db=None, n_threads=1, **doit_config_kwds):
        '''Stores tasks and the files they operate on, along with
        doit config and other metadata. This is the core of the pipelines:
        it passes its tasks along to doit for execution, and can check task
        and pipeline completion status.

        Args:
            directory (str): The directory in which to run the tasks. Will be
                created it it doesn't exist.
            logger (logging.Logger): Logger to record to.
            files (dict): Files used by the handler. Starts empty if omitted.
            profile (bool): If True, profile task execution.
            db (str): Name of the doit database.
            **doit_config_kwds: Keyword arguments passed to doit.

        Attributes:
            files (dict): Files used by the tasks.
            directory (str): Working directory for execution.
            tasks (OrderedDict): The tasks to execute.
            dep_file (str): Path of the doit database.
            doit_config (dict): The doit configuration given to the task runner.
            doit_dep_mgr (doit.dependency.Dependency): Doit object to track task
                status.
            profile (bool): Whether to run the profiler on tasks.
            logger (logging.Logger): Logger to use.
        '''

        super(TaskHandler, self).__init__()

        if files is None:
            self.files = {}
        elif type(files) is not dict:
            raise TypeError('files must be of type dict')
        else:
            self.files = files

        self.tasks = OrderedDict()
        
        self.directory = directory
        try:
            mkdir(directory)
        except OSError:
            pass

        if db is None:
            dep_file = path.join(self.directory, 'doit.db')
        else:
            dep_file = path.join(self.directory, '{0}.doit.db'.format(db))
        self.dep_file = dep_file
        logger.debug('Dependency Database File: {0}'.format(dep_file))
        self.doit_config = dict(dep_file=self.dep_file,
                                reporter=ui.GithubMarkdownReporter,
                                **doit_config_kwds)
        self.doit_dep_mgr = Dependency(SqliteDB, dep_file)
        self.n_threads = n_threads
        self.profile = profile
        self.logger = logger
        

    def register_task(self, name, task, files=None):
        '''Register a new task and its files with the handler.

        It may seem redundant or confusing to give the tasks a name different
        than their internal doit name. I do this because doit tasks need to have 
        names as unique as possible, so that they can be reused in different
        projects. A particular TaskHandler instance is only used for one
        pipeline run, and allowing different names makes it easier to reference
        tasks from elsewhere.

        Args:
            name (str): Name of the task. Does not have to correspond to doit's
                internal task name.
            task (:obj:): Either a dictionary or Task object.
            files (dict): Dictionary of files used.
        '''

        if files is None:
            files = {}
        if type(files) is not dict:
            raise TypeError('files must be of type dict')
        
        self.tasks[name] = task
        self.files.update(files)
        self.logger.debug('registered task {0}: {1}\n'
                          '  with files {2}'.format(name, task, files))

    def clear_tasks(self):
        '''Empty the task dictionary.'''

        self.logger.debug('Clearing {0} tasks'.format(len(self.tasks)))
        self.tasks = {}

    def get_status(self, task, move=False):
        '''Get the up-to-date status of a single task.

        Args:
            task (str): The task name to look up.
            move (bool): If True, move to the handler's directory before
                checking. Whether this is necessary depends mostly on whether
                the task uses relative or absolute paths.
        Returns:
            str: The string represenation of the status. Either "run" or
            "uptodate".
        '''

        if type(task) is str:
            try:
                task = self.tasks[task]
            except KeyError:
                self.logger.error('Task not found:{0}'.format(task))
                raise
        self.logger.debug('Getting status for task {0}'.format(task.name))
        if move:
            with Move(self.directory):
                status = self.doit_dep_mgr.get_status(task, self.tasks.values(),
                                                      get_log=True)
        else:
            status = self.doit_dep_mgr.get_status(task, self.tasks.values(),
                                                      get_log=True)
        self.logger.debug('Task {0} had status {1}'.format(task, status.status))
        try:
            self.logger.debug('Task {0} had reasons {1}'.format(task, status.reasons))
        except AttributeError:
            pass

        return status.status

    def print_statuses(self, uptodate_msg='All tasks up-to-date!',
                             outofdate_msg='Some tasks out of date!'):
        '''Print the up-to-date status of all tasks.

        Args:
            uptodate_msg (str): The message to print if all tasks are up to
            date.
        Returns:
            tuple: A bool (True if all up to date) and a dictionary of statuses.
        '''

        uptodate, statuses = self.check_uptodate()
        if uptodate:
            print(ui.paragraph(uptodate_msg))
        else:
            print(ui.paragraph(outofdate_msg))
            uptodate_list = [t for t,s in statuses.items() if s is True]
            outofdate_list = [t for t,s in statuses.items() if s is False]
            if uptodate_list:
                print('\nUp-to-date tasks:')
                print(ui.listing(uptodate_list))
            if outofdate_list:
                print('\nOut-of-date tasks:')
                print(ui.listing(outofdate_list))
        return uptodate, statuses

    def check_uptodate(self):
        '''Check if all tasks are up-to-date, ie if the pipeline is complete.
        Note that this moves to the handler's directory to lessen issues with
        relative versus absolute paths.

        Returns:
            bool: True if all are up to date.
        '''

        with Move(self.directory):
            statuses = {}
            outofdate = False
            for task_name, task in self.tasks.items():
                status = self.get_status(task)
                statuses[task_name] = status == 'up-to-date'
            return all(statuses.values()), statuses
        
    def load_tasks(self, cmd, opt_values, pos_args):
        '''Internal to doit -- triggered by the TaskLoader.'''

        self.logger.debug('loading {0} tasks'.format(len(self.tasks)))
        return self.tasks.values(), self.doit_config

    def run(self, doit_args=None, verbose=True):
        '''Run the pipeline. Movees to the directory, loads the tasks into doit,
        and executes that tasks that are not up-to-date.

        Args:
            doit_args (list): Args that would be passed to the doit shell
                command. By default, just run.
            verbose (bool): If True, print UI stuff.
        Returns:
            int: Exit status of the doit command.
        '''
        if verbose:
            print(ui.header('Run Tasks', level=4))
        if doit_args is None:
            doit_args = ['run']
            if self.n_threads > 1:
                doit_args.extend(['-n', str(self.n_threads)])

        runner = DoitMain(self)

        with Move(self.directory):
            if self.profile is True:
                profile_fn = path.join(self.directory, 'profile.csv')
                with StartProfiler(filename=profile_fn):
                    return runner.run(doit_args)
            else:
                return runner.run(doit_args)

