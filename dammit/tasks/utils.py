#!/usr/bin/env python
from __future__ import print_function

import os
from shutil import rmtree

from ..utils import doit_task

def clean_folder(target):
    try:
        rmtree(target)
    except OSError:
        pass


@doit_task
def get_group_task(group_name, tasks):

    return {'name': group_name,
            'actions': None,
            'task_dep': [t.name for t in tasks]}

