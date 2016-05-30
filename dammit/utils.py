#!/usr/bin/env python
from __future__ import print_function

import os
import sys

try:
    from pathlib import Path
except ImportError:
    from pathlib2 import Path

from doit.task import Task


class DammitTask(Task):

    def __repr__(self):
        return '{{ DammitTask: {name}'\
               '\n    actions: {actions}'\
               '\n   file_dep: {file_dep}'\
               '\n   task_dep: {task_dep}'\
               '\n    targets: {targets} }}'.format(actions=self.actions, **vars(self))


def dict_to_task(task_dict):
    if 'actions' not in task_dict:
        raise InvalidTask("Task %s must contain 'actions' field. %s" %
                          (task_dict['name'], task_dict))

    task_attrs = list(task_dict.keys())
    valid_attrs = set(Task.valid_attr.keys())
    for key in task_attrs:
        if key not in valid_attrs:
            raise InvalidTask("Task %s contains invalid field: '%s'"%
                              (task_dict['name'], key))

    return DammitTask(**task_dict)


def doit_task(task_dict_func):
    '''Wrapper to decorate functions returning pydoit
    Task dictionaries and have them return pydoit Task
    objects
    '''
    def d_to_t(*args, **kwargs):
        task_dict = task_dict_func(*args, **kwargs)
        return dict_to_task(task_dict)
    return d_to_t


class Move(object):

    def __init__(self, target, create=False):
        print('Move to', target, file=sys.stderr)
        self.target = target
        self.create = create
   
    def __enter__(self):
        self.cwd = os.getcwd()
        print('cwd:', self.cwd, file=sys.stderr)
        if self.create:
            try:
                os.mkdir(self.target)
            except OSError:
                pass
        os.chdir(self.target)

    def __exit__(self, exc_type, exc_value, traceback):
        os.chdir(self.cwd)
        if exc_type:
            return False


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


