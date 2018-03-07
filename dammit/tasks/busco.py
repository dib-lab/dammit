# Copyright (C) 2015-2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the BSD license.  See the LICENSE file for details.

import os

from doit.action import CmdAction
from doit.tools import run_once
from doit.task import clean_targets
import pandas as pd

from dammit import meta
from dammit.profile import profile_task
from dammit.utils import doit_task, which
from dammit.tasks.utils import clean_folder, DependentTask, InstallationError


class BuscoTask(DependentTask):

    def deps(self):
        buscov2 = which('BUSCO.py')
        buscov3 = which('run_BUSCO.py')

        tblastn = which('tblastn')
        makeblastdb = which('makeblastdb')
        if buscov2 is None and buscov3 is None:
            raise InstallationError('BUSCO not found. NOTE: '\
                                    'dammit 1.0 requires BUSCO v2 or greater')
        if tblastn is None:
            raise InstallationError('tblastn not found, required for BUSCO.')
        if makeblastdb is None:
            raise InstallationError('makeblastdb not found, required for BUSCO.')
        if self.logger:
            logger.debug('BUSCO:' + busco)
        return buscov3 if buscov3 is not None else buscov2

    @doit_task
    @profile_task
    def task(self, input_filename, output_name, busco_db_dir, 
                   input_type='tran', n_threads=1, params=None):
        '''Get a task to run BUSCO on the given FASTA file.

        Args:
            input_filename (str): The FASTA file to run BUSCO on.
            output_name (str): Base name for the BUSCO output directory.
            busco_db_dir (str): Directory with the BUSCO databases.
            input_type (str): By default, `trans` for transcriptome.
            n_threads (int): Number of threads to use.
            params (list): Extra parameters to pass to the executable.

        Returns:
            dict: A doit task.
        '''

        name = 'busco:{0}-{1}'.format(os.path.basename(input_filename),
                                      os.path.basename(busco_db_dir))

        exc = self.deps()
        config_file = os.path.join(meta.__path__, 'busco_config.ini')

        # BUSCO chokes on file paths as output names
        output_name = os.path.basename(output_name)
        cmd = ['BUSCO_CONFIG_FILE="{0}"'.format(config_file),
               'python3', exc, '-i', input_filename, '-f', '-o', output_name,
               '-l', busco_db_dir, '-m', input_type, '-c', str(n_threads)]
        if params is not None:
            cmd.extend(params)
        cmd = ' '.join(cmd)

        output_folder = 'run_' + output_name
        target_fn = os.path.join(output_folder, 'full_table_{0}.tsv'.format(output_name))

        return {'name': name,
                'actions': [cmd],
                'file_dep': [input_filename],
                'targets': [target_fn],
                'clean': [(clean_folder, [output_folder])]}


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
