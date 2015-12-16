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
                   get_maf_best_hits_task, \
                   get_annotate_fasta_task

def get_report_tasks(transcriptome, annotator, databases, n_threads=1):

    tasks = []
    outputs = []


    orthodb_best_hits = annotator.orthodb_fn + '.best.csv'
    orthodb_gff3 = annotator.orthodb_fn + '.gff3'
    tasks.append(
        get_maf_best_hits_task(annotator.orthodb_fn,
                               orthodb_best_hits)
    )
    tasks.append(
        get_maf_gff3_task(orthodb_best_hits,
                          orthodb_gff3, 'OrthoDB')
    )
    outputs.append(orthodb_gff3)

    pfam_gff3 = annotator.pfam_fn + '.gff3'
    tasks.append(
        get_hmmscan_gff3_task(annotator.pfam_fn,
                              pfam_gff3, 
                              'Pfam')
    )
    outputs.append(pfam_gff3)

    rfam_gff3 = annotator.rfam_fn + '.gff3'
    tasks.append(
        get_cmscan_gff3_task(annotator.rfam_fn,
                             rfam_gff3, 
                             'Rfam')
    )
    outputs.append(rfam_gff3)

    for db, fn in annotator.user_pep_fn_dict.iteritems():
        gff3_fn = fn + '.gff3'
        tasks.append(
            get_crb_gff3_task(fn, gff3_fn, db)
        )
        outputs.append(gff3_fn)

    outputs.append(annotator.transdecoder_gff3_fn)

    merged_gff3 = transcriptome + '.dammit.gff3'
    tasks.append(
        get_gff3_merge_task(outputs, merged_gff3)
    )
    outputs.append(merged_gff3)

    renamed_transcriptome = transcriptome + '.dammit.fasta'
    tasks.append(
        get_annotate_fasta_task(transcriptome, merged_gff3,
                                renamed_transcriptome)
    )
    outputs.append(renamed_transcriptome)

    return outputs, tasks

