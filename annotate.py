#!/usr/bin/env python
from __future__ import print_function

import os
from platform import system
import sys

import common
from tasks import get_blast_format_task, \
                  get_transcriptome_stats_task, \
                  get_busco_task, \
                  get_group_task, \
                  get_link_file_task

def get_annotate_tasks(transcriptome, prog_paths, database_dict, 
                       n_threads=1, user_databases=[]):

    tasks = []

    user_database_dict = {}
    if user_databases:
        for db in user_databases:
            db_path = db + '.db'
            tasks.append(
                get_blast_format_task(db, db_path, 'prot')
            )
            user_databases[os.path.basename(db)] = db_path

    tasks.append(
            get_link_file_task(os.path.abspath(transcriptome))
    )

    '''
    Calculate assembly information. First it runs some basic stats like N50 and
    number of contigs, and uses the HyperLogLog counter from khmer to
    estimate unique k-mers for checking redundancy. Then it runs BUSCO to
    assess completeness. These tasks are grouped under the 'assess' task.
    '''
    assess_tasks = []
    assess_tasks.append(
        get_transcriptome_stats_task(transcriptome, 
                                     os.path.basename(transcriptome + '.stats'))
    )

    busco_cfg = common.CONFIG['settings']['busco']
    busco_output_name = '{0}.busco.results'.format(transcriptome)
    assess_tasks.append(
        get_busco_task(transcriptome, busco_output_name, database_dict['BUSCO'],
                       'trans', n_threads, busco_cfg,
                       busco_dir=prog_paths['busco'])
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

    cwd = os.getcwd()
    try:
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        os.chdir(output_dir)

        common.run_tasks(tasks, args, config=doit_config)
    finally:
        os.chdir(cwd)

