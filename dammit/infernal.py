#!/usr/bin/env python
import os

from doit.action import CmdAction
from doit.tools import title_with_actions
from doit.task import clean_targets

from .common import which
from .parallel import parallel_fasta, multinode_parallel_fasta
from .tasks import doit_task


@doit_task
def get_cmpress_task(db_filename, infernal_cfg):

    exc = which('cmpress')
    cmd = '{exc} {db_filename}'.format(**locals())

    return {'name': 'cmpress:' + os.path.basename(db_filename),
            'title': title_with_actions,
            'actions': [cmd],
            'targets': [db_filename + ext for ext in ['.i1f', '.i1i', '.i1m', '.i1p']],
            'uptodate': [True],
            'clean': [clean_targets]}


@doit_task
def get_cmscan_task(input_filename, output_filename, db_filename,
                    cutoff, n_threads, infernal_cfg, n_nodes=None):

    name = 'cmscan:' + os.path.basename(input_filename) + '.x.' + \
           os.path.basename(db_filename)

    exc = which('cmscan')
    if n_nodes is None:
        parallel_cmd = parallel_fasta(input_filename, n_threads)
    else:
        parallel_cmd = multinode_parallel_fasta(input_filename, n_threads,
                                                n_nodes)

    stat = output_filename + '.out'
    cmd = [parallel_cmd, exc, '--cpu', '1', '--rfam', '--nohmmonly',
           '-E', str(cutoff), '--tblout', '/dev/stdout', '-o', stat,
           db_filename, '/dev/stdin', '>', output_filename]
    cmd = ' '.join(cmd)

    return {'name': name,
            'title': title_with_actions,
            'actions': [cmd],
            'file_dep': [input_filename, db_filename, db_filename + '.i1p'],
            'targets': [output_filename, output_filename + '.cmscan'],
            'clean': [clean_targets]}


