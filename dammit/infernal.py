#!/usr/bin/env python
import os

from doit.action import CmdAction
from doit.tools import title_with_actions, LongRunning
from doit.task import clean_targets

from .utils import which, doit_task, convert_pathlib
from .parallel import parallel_fasta, multinode_parallel_fasta, get_filesize_task


@doit_task
@convert_pathlib
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
@convert_pathlib
def get_cmscan_task(input_filename, output_filename, db_filename,
                    cutoff=0.00001, n_threads=1, n_nodes=None, params=None):

    name = 'cmscan:' + os.path.basename(input_filename) + '.x.' + \
           os.path.basename(db_filename)
    stat = output_filename + '.cmscan.out'

    def cmscan_cmd(file_size):
        exc = which('cmscan')
        if n_nodes is None:
            parallel_cmd = parallel_fasta(input_filename, n_threads, file_size)
        else:
            parallel_cmd = multinode_parallel_fasta(input_filename, n_threads,
                                                    n_nodes, file_size)


        cmd = [parallel_cmd, exc]
        if params is not None:
            cmd.extend([str(p) for p in params])
        cmd.extned(['--cpu', '1', '--rfam', '--nohmmonly',
               '-E', str(cutoff), '--tblout', '/dev/stdout', '-o', stat,
               db_filename, '/dev/stdin', '>', output_filename])
        return ' '.join(cmd)

    return {'name': name,
            'title': title_with_actions,
            'actions': [LongRunning(cmscan_cmd)],
            'file_dep': [input_filename, db_filename, db_filename + '.i1p'],
            'targets': [output_filename, output_filename + '.cmscan'],
            'getargs': {'file_size': ('get_filesize:{0}'.format(input_filename),
                                      'size')},
            'clean': [clean_targets]}


def cmscan(input_filename, output_filename, db_filename, **cmscan_kwds):



    yield get_filesize_task(input_filename)
    yield get_cmscan_task(input_filename, output_filename, db_filename,
                          **cmscan_kwds)
