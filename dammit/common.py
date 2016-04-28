#!/usr/bin/env python
from __future__ import print_function

import json
import logging
import os
import sys
import textwrap
import time

from doit.cmd_base import TaskLoader
from doit.doit_cmd import DoitMain

CUR_TIME = time.strftime('%Y-%m-%d-%H%M')

def uniqify_tasks(tasks):
    seen = set()
    uniqified = []
    for task in tasks:
        if task.name not in seen:
            uniqified.append(task)
            seen.add(task.name)
    return uniqified


