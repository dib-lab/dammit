# Copyright (C) 2015-2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the BSD license.  See the LICENSE file for details.

from functools import wraps
import os
import stat
import sys

from doit.action import PythonAction
from doit.task import Task, InvalidTask
import six


def cleaned_actions(actions):
    '''Get a cleanup list of actions: Python actions
    have their <locals> portion stripped, which clutters 
    up PythonActions that are closures.
    '''
    txt = ''
    for action in actions:
        txt_rep = six.text_type(action)
        if isinstance(action, PythonAction):
            # clean up inner fuctions in Python actions
            txt_rep = txt_rep.replace('<locals>.', '')
        else:
            txt_rep = txt_rep[:5] + '`' + txt_rep[5:] + '`'
        txt += "\n    * {0}".format(txt_rep)
    return txt


class DammitTask(Task):
    '''Subclass doit.task.Task for dammit. Updates the string __repr__
    and adds a uniform updated title function.
    '''

    def __repr__(self):
        return '{{ DammitTask: {name}'\
               '\n    actions: {actions}'\
               '\n   file_dep: {file_dep}'\
               '\n   task_dep: {task_dep}'\
               '\n    targets: {targets} }}'.format(actions=self.actions,
                                                    **vars(self))

    def title(self):
        if self.custom_title:
            return self.custom_title(self)
        else:
            if self.actions:
                title = cleaned_actions(self.actions)
            else:
                title = "Group: %s" % ", ".join(self.task_dep)
            return "%s: %s"% (self.name, title)


def dict_to_task(task_dict):
    '''Given a doit task dict, return a DammitTask.

    Args:
        task_dict (dict): A doit task dict.

    Returns:
        DammitTask: Subclassed doit task.
    '''

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
    @wraps(task_dict_func)
    def d_to_t(*args, **kwargs):
        task_dict = task_dict_func(*args, **kwargs)
        return dict_to_task(task_dict)
    return d_to_t


def touch(filename):
    '''Perform the equivalent of bash's touch on the file.

    Args:
        filename (str): File path to touch.
    '''

    open(filename, 'a').close()


class Move(object):
    '''Context manager to change current working directory.
    '''

    def __init__(self, target, create=False, verbose=False):
        '''Move to specified directory.

        Args:
            target (str): Directory to change to.
            create (bool): If True, create the directory.
        '''

        self.verbose = verbose
        self.target = target
        self.create = create
   
    def __enter__(self):
        self.cwd = os.getcwd()
        if self.verbose:
            print('Move to `{0}` from cwd: `{1}`'.format(self.target, 
                                                     self.cwd, 
                                                     file=sys.stderr))
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


