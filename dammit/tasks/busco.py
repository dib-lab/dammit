#!/usr/bin/env python
from __future__ import print_function

import os

from doit.action import CmdAction
from doit.tools import run_once, title_with_actions
from doit.task import clean_targets
import pandas as pd

from .utils import clean_folder
from ..profile import profile_task
from ..utils import doit_task, which


@doit_task
@profile_task
def get_busco_task(input_filename, output_name, busco_db_dir, 
                   input_type='trans', n_threads=1, params=None):

    name = 'busco:{0}-{1}'.format(os.path.basename(input_filename),
                                  os.path.basename(busco_db_dir))

    assert input_type in ['genome', 'OGS', 'trans']
    exc = which('BUSCO_v1.1b1.py')
    # BUSCO chokes on file paths as output names
    output_name = os.path.basename(output_name)
    cmd = ['python3', exc, '-in', input_filename, '-f', '-o', output_name,
           '-l', busco_db_dir, '-m', input_type, '-c', str(n_threads)]
    if params is not None:
        cmd.extend(params)
    cmd = ' '.join(cmd)

    return {'name': name,
            'title': title_with_actions,
            'actions': [cmd],
            'file_dep': [input_filename],
            'uptodate': [run_once],
            'clean': [(clean_folder, ['run_' + output_name])]}


def parse_busco_full(fn):
    '''Parses a BUSCO full result table into a Pandas DataFrame.

    Args:
        fn (str): The results file.
    Returns:
        DataFrame: The results DataFrame.
    '''

    df = pd.read_table(fn)
    return df.rename(columns={'#BUSCO_group': 'BUSCO_group'})


def parse_busco_summary(fn):
    '''Parses a BUSCO summary file into a JSON compatible
    dictionary.

    Args:
        fn (str): The summary results file.
    Returns:
        dict: The BUSCO results.
    '''

    res = {}
    with open(fn) as fp:
        for ln in fp:
            if ln.strip().startswith('C:'):
                tokens = ln.split(',')
                for token in tokens:
                    key, _, val = token.partition(':')
                    key = key.strip()
                    val = val.strip().strip('%')
                    if key == 'C':
                        valc, _, vald = val.partition('%')
                        valc = valc.strip()
                        vald = vald.strip('D:][%')
                        res['C(%)'] = valc
                        res['D(%)'] = vald
                    else:
                        if key != 'n':
                           key += '(%)'
                        res[key] = val.strip().strip('%')
    return res


def parse_busco_multiple(fn_list, dbs=['metazoa', 'vertebrata']):
    '''Parses multiple BUSCO results summaries into an appropriately
    index DataFrame.

    Args:
        fn_list (list): List of paths to results files.
        dbs (list): List of BUSCO database names.
    Returns:
        DataFrame: The formated DataFrame.
    '''

    data = []
    for fn in fn_list:
        data.append(parse_busco_summary(fn))

    df = pd.DataFrame(data)
    df['fn'] = [os.path.basename(fn)[14:-14].strip('.') for fn in fn_list]
    df['db'] = None
    for db in dbs:
        idx = df.fn.str.contains(db)
        df.loc[idx,'db'] = db
        df.loc[idx,'fn'] = df.loc[idx, 'fn'].apply(lambda fn: fn[:fn.find(db)].strip('. '))

    return df


def busco_to_df(fn_list, dbs=['metazoa', 'vertebrata']):
    ''' Given a list of BUSCO results from different databases, produce
    an appropriately multi-indexed DataFrame of the results.

    Args:
        fn_list (list): The BUSCO summary files.
        dbs (list): The BUSCO databases used for these runs.
    Returns:
        DataFrame: The BUSCO results.
    '''

    data = []
    for fn in fn_list:
        data.append(parse_busco(fn))

    df = pd.DataFrame(data)
    df['fn'] = [os.path.basename(fn)[14:-14].strip('.') for fn in fn_list]
    df['db'] = None
    for db in dbs:
        idx = df.fn.str.contains(db)
        df.loc[idx,'db'] = db
        df.loc[idx,'fn'] = df.loc[idx, 'fn'].apply(lambda fn: fn[:fn.find(db)].strip('. '))
    return df
