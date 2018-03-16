# Copyright (C) 2015-2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the BSD license.  See the LICENSE file for details.

import os
from shutil import rmtree

from dammit.utils import doit_task

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

