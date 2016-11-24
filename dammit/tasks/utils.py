#!/usr/bin/env python
from __future__ import print_function

import os
from shutil import rmtree

from ..utils import doit_task

def clean_folder(target):
    '''Function for doit task's `clean` parameter to remove a folder.

    Args:
        target (str): The folder to remove.
    '''

    try:
        rmtree(target)
    except OSError:
        pass

class DependentTask:

    def __init__(self, logger=None):
        self.logger = logger

    def deps(self):
        raise NotImplementedError()

    def task(self, *args, **kwargs):
        raise NotImplementedError()


class InstallationError(RuntimeError):
    pass


@doit_task
def get_group_task(group_name, tasks):
    '''Creat a task group from the given tasks.

    Args:
        group_name (str): The name to give the group.
        tasks (list): List of Task objects to add to group.

    Returns:
        dict: A doit task for the group.
    '''

    return {'name': group_name,
            'actions': None,
            'task_dep': [t.name for t in tasks]}

