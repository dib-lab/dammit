#!/usr/bin/env python
from __future__ import print_function

import os
import sys

from . import common
from .tasks import get_gff3_report_task

def get_report_tasks(transcriptome, results_dict, n_threads=1):

    tasks = []
    results = {}

    tasks.append(
        get_gff3_report_task(transcriptome, results_dict)
    )

    return tasks

