# Copyright (C) 2015-2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the BSD license.  See the LICENSE file for details.

import csv
import filelock
from functools import wraps
import sys
import time
import warnings
from contextlib import contextmanager
from os import path

from doit.task import Task as DoitTask

from dammit.utils import cleaned_actions


class Profiler(object):
    '''Thread-safe performance profiler.
    '''

    def __init__(self):
        self.running = False

    def start_profiler(self, filename=None, blockname='__main__'):
        '''Start the profiler, with results stored in the given filename.

        Args:
            filename (str): Path to store profiling results. If not given,
                uses a representation of the current time
            blockname (str): Name assigned to the main block.
        '''

        self.run_name = time.strftime("%a_%d_%b_%Y_%H%M%S", time.localtime())
        if filename is None:
            self.filename = '{0}.csv'.format(self.run_name)
        else:
            self.filename = filename
        self.run_name = time.ctime()
        self.start_time = time.time()
        self.blockname = blockname
        self.running = True
        self.lock = filelock.FileLock('{0}.lock'.format(self.filename))
        print('Profiling is ON:', self.filename, '\n', file=sys.stderr)

    def write_result(self, task_name, start_time, end_time, elapsed_time):
        '''Write results to the file, using the given task name as the
        name for the results block.

        Args:
            task_name (str): ID for the result row (the block profiled).
            start_time (float): Time of block start.
            end_time (float): Time of block end.
            elapsed_time (float): Total time.
        '''

        try:
            with self.lock.acquire(timeout=10):
                header = not path.isfile(self.filename)
                with open(self.filename, 'a') as fp:
                    writer = csv.writer(fp, delimiter=',')
                    if header:
                        writer.writerow(['run_id', 'block', 'start_t', 'end_t',
                                      'elapsed_t'])
                    row = [self.run_name, task_name, start_time, end_time, elapsed_time]
                    writer.writerow(row)
        except filelock.Timeout as e:
            warnings.warn(e, RuntimeWarning, stacklevel=1)

    def stop_profiler(self):
        '''Shut down the profiler and write the final elapsed time.
        '''
        self.end_time = time.time()
        elapsed = self.end_time - self.start_time
        self.write_result(self.blockname, self.start_time, self.end_time, elapsed)
        self.running = False
        return elapsed


class Timer(object):
    '''Simple timer class.
    '''

    def start(self):
        '''Start the timer.
        '''
        self.start_time = time.time()

    def stop(self):
        '''Stop the timer and return the elapsed time.
        '''
        self.end_time = time.time()
        return self.end_time - self.start_time


def title_without_profile_actions(task):
    """Generate title without profiling actions"""
    title = ''
    if task.actions:
        title = cleaned_actions(task.actions[1:-1])
    else:
        title = "Group: %s" % ", ".join(task.task_dep)
    return "%s: %s"% (task.name, title)


def setup_profiler():
    '''Returns a context manager, a funnction to add profiling actions to doit
    tasks, and a decoratator to apply that function to task functions.

    The profiling function adds new actions to the beginning and end of
    the given task's action list, which start and stop the profiler and 
    record the results. The task decorator applies this function. The actions
    only record data if the profiler is running when they are called, and
    they are removed from doit's execution output to reduce clutter.

    The context manager starts the profiler in its block, storing data
    in the given file.

    Yes, this is a function function function which creates six different
    functions at seven different function scopes. Written in honor
    of javascript programmers everywhere, and to baffle and irritate
    @ryneches.
    '''

    profiler = Profiler()
    @contextmanager
    def profiler_manager(filename=None, blockname='__main__'):
        profiler.start_profiler(filename=filename, blockname=blockname)
        yield
        profiler.stop_profiler()

    def add_profile_actions(task):
        timer = Timer()
        def start_profiling():
            if profiler.running:
                timer.start()

        def stop_profiling():
            if profiler.running:
                elapsed = timer.stop()
                profiler.write_result(task['name'], timer.start_time,
                                      timer.end_time, elapsed)
        
        if isinstance(task, DoitTask):
            actions = task._actions
            task.custom_title = title_without_profile_actions
        else:
            actions = task['actions']
            task['title'] = title_without_profile_actions

        actions.insert(0, start_profiling)
        actions.append(stop_profiling)

        return task
    
    def profile_decorator(task_func):

        @wraps(task_func)
        def func(*args, **kwargs):
            task = task_func(*args, **kwargs)
            return add_profile_actions(task)
        
        return func

    return profiler_manager, add_profile_actions, profile_decorator

StartProfiler, add_profile_actions, profile_task = setup_profiler()

