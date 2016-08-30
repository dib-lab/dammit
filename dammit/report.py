#!/usr/bin/env python
from __future__ import print_function

import os
import sys


def get_report_tasks(transcriptome, annotator, databases, n_threads=1):

    tasks = []
    outputs = []


    for db, fn in annotator.user_pep_fn_dict.iteritems():
        gff3_fn = fn + '.gff3'
        tasks.append(
            get_crb_gff3_task(fn, gff3_fn, db)
        )
        outputs.append(gff3_fn)

    return outputs, tasks

