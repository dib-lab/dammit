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
