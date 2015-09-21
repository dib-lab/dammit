#!/usr/bin/env python
from __future__ import print_function

import os
from platform import system
import sys

import common
from tasks import get_blast_format_task

def get_annotate_tasks(transcriptome, output_dir, prog_paths, databases):

    tasks = []

    user_databases = []
    if args.user_databases:
        for db in args.databases:
            db_name = os.path.join(out_dir, db + '.db')
            tasks.append(
                get_blast_format_task(db, db_name, 'prot')
            )
            user_databases.append(db_name)

    '''
    Calculate assembly information. First it runs some basic stats like N50 and
    number of contigs, and uses the HyperLogLog counter from khmer to
    estimate unique k-mers for checking redundancy. Then it runs BUSCO to
    assess completeness. These tasks are grouped under the 'assess' task.
    '''
    assess_tasks = []
    assess_tasks.append(
        get_transcriptome_stats_task(transcriptome, out_dir)
    )

    busco_cfg = common.CONFIG['settings']['busco']
    busco_dir = os.path.join(output_dir,
                             '{0}.busco.results'.format(transcriptome))
    assess_tasks.append(
        busco_task(transcriptome, busco_dir, databases['BUSCO'], 
                   'trans', busco_cfg)
    )

    tasks.extend(assess_tasks)
    tasks.append(get_group_task('assess', assess_tasks))


    return tasks

def run_annotate_tasks(transcriptome, output_dir, tasks, args=['run']):
    '''
    Set up doit's config for the actual analysis tasks.
    We'll put the doit database for these tasks into the output
    directory so that we don't end up scattering them around the
    filesystem, or worse, with one master db containing dependency
    metadata from every analysis ever run by the user!
    '''

    doit_config = {
                    'backend': common.DOIT_BACKEND,
                    'verbosity': common.DOIT_VERBOSITY,
                    'dep_file': os.path.join(output_dir, '.' +
                                             os.path.basename(transcriptome) +
                                             '.doit.db')
                  }
    run_tasks(tasks, args, config=doit_config)
