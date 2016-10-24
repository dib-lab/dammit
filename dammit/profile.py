import csv
import filelock
import sys
import time
import warnings
from contextlib import contextmanager
from os import path
import six

from doit.task import Task as DoitTask


class Profiler(object):

    def __init__(self):
        self.running = False

    def start_profiler(self, filename=None, blockname='__main__'):
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
        self.end_time = time.time()
        elapsed = self.end_time - self.start_time
        self.write_result(self.blockname, self.start_time, self.end_time, elapsed)
        self.running = False
        return elapsed


class Timer(object):

    def start(self):
        self.start_time = time.time()

    def stop(self):
        self.end_time = time.time()
        return self.end_time - self.start_time


def title_without_profile_actions(task):
    """return task name task actions"""
    title = ''
    if task.actions:
        for action in task.actions[1:-1]:
            title += "\n\t* " + six.text_type(action)
    # A task that contains no actions at all
    # is used as group task
    else:
        title = "Group: %s" % ", ".join(task.task_dep)
    return "%s: %s"% (task.name, title)


def setup_profiler():
    '''Returns a doit task decorator and a context manager for profiling.

    The task decorator calls records profile data when the task finishes
    executing by adding new tasks to the beginning and end of the actions
    list. The context manager simply starts a new profiler instance within
    its block.

    Yes, this is a function function funcion which creates four different
    functions at seven different function scopes. Written in honor
    of javascript programmers everywhere.
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

        def func(*args, **kwargs):
            task = task_func(*args, **kwargs)
            return add_profile_actions(task)
        
        return func

    return profiler_manager, add_profile_actions, profile_decorator

StartProfiler, add_profile_actions, profile_task = setup_profiler()

