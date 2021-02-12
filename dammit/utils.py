# Copyright (C) 2015-2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the BSD license.  See the LICENSE file for details.

import os
import sys
import yaml
import collections.abc

import click
import numpy as np
from ope.io.base import ChunkParser
import pandas as pd


class ShortChoice(click.Choice):
    ''' Modified click.Choice parameter type that truncates
    the list of choices.
    '''

    def get_metavar(self, param):
        return f"[{'|'.join(self.choices[:5])}|...]"


def update_nested_dict(d, other):
    # Can't just update at top level, need to update nested params
    # Note that this only keeps keys that already exist in other
    # https://code.i-harness.com/en/q/3154af
    for k, v in other.items():
        if isinstance(v, collections.abc.Mapping):
            d_v = d.get(k)
            if isinstance(d_v, collections.abc.Mapping):
                update_nested_dict(d_v, v)
            else:
                d[k] = v.copy()
        else:
            d[k] = v


def create_dirs(dirs):
    for d in dirs:
        os.makedirs(d, exist_ok=True)


def touch(filename):
    '''Perform the equivalent of bash's touch on the file.

    Args:
        filename (str): File path to touch.
    '''

    open(filename, 'a').close()


class Namespace:
    pass


def read_yaml(filename):
    with open(filename, 'r') as stream:
        try:
            yamlD = yaml.safe_load(stream) #, Loader=yaml.FullLoader)
        except yaml.YAMLError as exc:
            print(exc, file=sys.stderr)
            print(f'ERROR parsing YAML from {filename}, exiting.', file=sys.stderr)
            sys.exit(1)
    return yamlD


def write_yaml(yamlD, paramsfile):
    with open(paramsfile, 'w') as params:
        yaml.dump(yamlD, stream=params, indent=2, default_flow_style=False)


class BUSCOTableParser(ChunkParser):
    
    columns = [('BUSCO_id', str),
               ('Status', str),
               ('Sequence', str),
               ('Score', float),
               ('Length', int),
               ('Start', int),
               ('End', int)]
    
    def __init__(self, filename, **kwargs):
        super(BUSCOTableParser, self).__init__(filename, **kwargs)
    
    def __iter__(self):
        with open(self.filename) as fp:
            header = fp.readline()
        self.busco_version = header.partition(':')[-1].strip()

        df = pd.read_table(self.filename,
                           names=[k for k,_ in self.columns[:-2]],
                           delimiter='\t',
                           comment='#',
                           error_bad_lines=False)
        
        if self.busco_version == '5.0.0' and not df.Sequence.isna().all():
            seq_df = df.Sequence.str.partition(':')
            coords_df = seq_df[2].str.partition('-')
            df['Sequence'] = seq_df[0]
            df['Start'] = pd.to_numeric(coords_df[0])
            df['End'] = pd.to_numeric(coords_df[2])
        else:
            df['Start'] = np.nan
            df['End'] = np.nan
        
        setattr(df, 'busco_version', self.busco_version)

        yield df.convert_dtypes()
