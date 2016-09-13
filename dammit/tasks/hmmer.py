#!/usr/bin/env python
from __future__ import print_function

from doit.action import CmdAction
from doit.tools import title_with_actions
from doit.task import clean_targets
import os
import pandas as pd

from .parallel import parallel_fasta
from . import parsers
from .utils import doit_task, which


@doit_task
def get_hmmscan_task(input_filename, output_filename, db_filename,
                     cutoff=0.00001, n_threads=1, pbs=False, 
                     params=None):

    name = 'hmmscan:' + os.path.basename(input_filename) + '.x.' + \
                    os.path.basename(db_filename)
    stat = output_filename + '.hmmscan.out'
    
    hmmscan_exc = which('hmmscan')
    cmd = [hmmscan_exc]
    if params is not None:
        cmd.extend([str(p) for p in params])
    cmd.extend(['--cpu', '1', '--domtblout', '/dev/stdout', 
                '-E', str(cutoff), '-o', stat, db_filename, '/dev/stdin'])
    
    cmd = parallel_fasta(input_filename, output_filename, cmd, n_threads, pbs=pbs)
        
    return {'name': name,
            'title': title_with_actions,
            'actions': [cmd],
            'file_dep': [input_filename, db_filename, db_filename+'.h3p'],
            'targets': [output_filename, stat],
            'clean': [clean_targets]}


@doit_task
def get_hmmpress_task(db_filename, params=None, task_dep=None):

    name = 'hmmpress:' + os.path.basename(db_filename)
    exc = which('hmmpress')

    cmd = [exc]
    if params is not None:
        cmd.extend([str(p) for p in params])
    cmd.append(db_filename)

    cmd = ' '.join(cmd)

    task_d =  {'name': name,
              'title': title_with_actions,
              'actions': [cmd],
              'targets': [db_filename + ext for ext in ['.h3f', '.h3i', '.h3m', '.h3p']],
              'uptodate': [True],
              'clean': [clean_targets]}
    if task_dep is not None:
        task_d['task_dep'] = task_dep

    return task_d


@doit_task
def get_remap_hmmer_task(hmmer_filename, remap_gff_filename, output_filename):

    name = 'remap_hmmer:{0}'.format(os.path.basename(hmmer_filename))

    def cmd():
        gff_df = pd.concat(parsers.parse_gff3(remap_gff_filename))
        hmmer_df = pd.concat(parsers.hmmscan_to_df_iter(hmmer_filename))

        merged_df = pd.merge(hmmer_df, gff_df, left_on='full_query_name', right_on='ID')

        hmmer_df['env_coord_from'] = (merged_df.start + \
                                      (3 * merged_df.env_coord_from)).astype(int)
        hmmer_df['env_coord_to'] = (merged_df.start + \
                                    (3 * merged_df.env_coord_to)).astype(int)
        hmmer_df['ali_coord_from'] = (merged_df.start + \
                                      (3 * merged_df.ali_coord_from)).astype(int)
        hmmer_df['ali_coord_to'] = (merged_df.start + \
                                    (3 * merged_df.ali_coord_to)).astype(int)

        hmmer_df.to_csv(output_filename, header=True, index=False)

    return {'name': name,
            'title': title_with_actions,
            'actions': [cmd],
            'file_dep': [hmmer_filename, remap_gff_filename],
            'targets': [output_filename],
            'clean': [clean_targets]}

