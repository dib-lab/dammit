#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2020
# File   : config.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 23.04.2020

import os
import shutil
import sys
import yaml

import click

from ..config import USER_CONFIG


@click.group(name='config')
@click.pass_obj
def config_group(config):
    ''' Show dammit configuration information.'''
    pass


@config_group.command('show-directories')
@click.pass_obj
@click.option('--database-dir',
              envvar='DAMMIT_DB_DIR',
              hidden=True)
@click.option('--temp-dir',
              envvar='DAMMIT_TEMP_DIR',
              hidden=True)
@click.option('--conda-dir',
              envvar='DAMMIT_CONDA_DIR',
              help='Directory to store snakemake-created conda environments.')
@click.option('--save', default='-', type=click.File('w'))
def show_directories_cmd(config, database_dir, temp_dir, conda_dir, save):
    '''List dammit directory locations.
    
    Locations come either from defaults or environment variables.'''

    if database_dir:
        config.core['database_dir'] = database_dir

    if temp_dir:
        config.core['temp_dir'] = temp_dir
    
    if conda_dir:
        config.core['conda_env_dir'] = conda_dir

    print('Databases dir:', config.core['database_dir'], file=save)
    print('Temp dir:', config.core['temp_dir'], file=save)
    print('Conda env dir:', config.core['conda_env_dir'], file=save)


@config_group.command('show-default')
@click.pass_obj
@click.argument('config', type=click.Choice(['user', 'core', 'databases', 'pipelines']))
@click.option('--save', default='-', type=click.File('w'))
def show_default_cmd(defaults, config, save):
    ''' Show the selected default configuration file.'''

    if config == 'user':
        save.write(yaml.dump(USER_CONFIG))
    else:
        save.write(yaml.dump(getattr(defaults, config)))


@config_group.command('busco-groups')
@click.pass_obj
@click.option('--save', default='-', type=click.File('w'))
def busco_groups_cmd(defaults, save):
    ''' Lists the available BUSCO group databases.'''

    save.write(' '.join(defaults.databases['busco']['lineages']))
