#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2020
# File   : run.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 23.04.2020

import os

import pandas as pd


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
