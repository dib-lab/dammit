#!/usr/bin/env python
import os

from doit.action import CmdAction
from doit.tools import title_with_actions
from doit.task import clean_targets

from .utils import which, doit_task
from .parallel import parallel_fasta


@doit_task
def get_cmpress_task(db_filename, params=None, task_dep=None):

    cmd = [which('cmpress')]
    if params is not None:
        cmd.extend([str(p) for p in params])
    cmd.append(db_filename)
    cmd = ' '.join(cmd)

    task_d = {'name': 'cmpress:' + os.path.basename(db_filename),
              'title': title_with_actions,
              'actions': [cmd],
              'targets': [db_filename + ext for ext in ['.i1f', '.i1i', '.i1m', '.i1p']],
              'uptodate': [True],
              'clean': [clean_targets]}
    if task_dep is not None:
        task_d['task_dep'] = task_dep

    return task_d


@doit_task
def get_cmscan_task(input_filename, output_filename, db_filename,
                    cutoff=0.00001, n_threads=1, pbs=False, params=None):

    name = 'cmscan:' + os.path.basename(input_filename) + '.x.' + \
           os.path.basename(db_filename)
    stat = output_filename + '.cmscan.out'

    exc = which('cmscan')
    cmd = [exc]
    if params is not None:
        cmd.extend([str(p) for p in params])
    cmd.extend(['--cpu', '1', '--rfam', '--nohmmonly',
           '-E', str(cutoff), '--tblout', '/dev/stdout', '-o', stat,
           db_filename, '/dev/stdin'])
    cmd = parallel_fasta(input_filename, output_filename, cmd, n_threads,
                         pbs=pbs)

    return {'name': name,
            'title': title_with_actions,
            'actions': [cmd],
            'file_dep': [input_filename, db_filename, db_filename + '.i1p'],
            'targets': [output_filename, stat],
            'clean': [clean_targets]}
