#!/usr/bin/env python
from __future__ import print_function

import os
import sys

from . import common
from .tasks import get_maf_gff3_task, \
                   get_hmmscan_gff3_task, \
                   get_cmscan_gff3_task, \
                   get_gff3_merge_task, \
                   get_crb_gff3_task, \
                   get_maf_best_hits_task

def get_report_tasks(transcriptome, results_dict, databases, n_threads=1):

    tasks = []
    outputs = []


    orthodb_best_hits = results_dict['orthodb'] + '.best.csv'
    orthodb_gff3 = results_dict['orthodb'] + '.gff3'
    tasks.append(
        get_maf_best_hits_task(results_dict['orthodb'],
                               orthodb_best_hits)
    )
    tasks.append(
        get_maf_gff3_task(orthodb_best_hits,
                          orthodb_gff3, 'OrthoDB')
    )
    outputs.append(orthodb_gff3)

    pfam_gff3 = results_dict['pfam'] + '.gff3'
    tasks.append(
        get_hmmscan_gff3_task(results_dict['pfam'],
                              pfam_gff3, 'Pfam')
    )
    outputs.append(pfam_gff3)

    rfam_gff3 = results_dict['rfam'] + '.gff3'
    tasks.append(
        get_cmscan_gff3_task(results_dict['rfam'],
                             rfam_gff3, 'Rfam')
    )
    outputs.append(rfam_gff3)

    for db, fn in results_dict['user'].iteritems():
        gff3_fn = fn + '.gff3'
        tasks.append(
            get_crb_gff3_task(fn, gff3_fn, db)
        )
        outputs.append(gff3_fn)

    outputs.append(results_dict['prot_predictions_gff3'])

    merged_gff3 = transcriptome + '.dammit.gff3'
    tasks.append(
        get_gff3_merge_task(outputs, merged_gff3)
    )
    outputs.append(merged_gff3)

    return outputs, tasks

