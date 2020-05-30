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

from ..config import (parse_config, DEFAULT_DATABASES_DIR, CONDA_ENV_TEMPDIR,
                      DEFAULT_TEMP_DIR, TEMP_SUBDIRS, WORKFLOW_CONFIG_TEMPDIR)


@click.group(name='config')
@click.pass_context
def config_group(ctx):
    ''' Show dammit configuration information.'''

    core, databases, pipelines = parse_config()
    ctx.obj = {'core': core, 'databases': databases, 'pipelines': pipelines}


@config_group.command('show-directories')
@click.option('--database-dir',
              default=DEFAULT_DATABASES_DIR,
              envvar='DAMMIT_DB_DIR',
              hidden=True)
@click.option('--temp-dir',
              default=DEFAULT_TEMP_DIR,
              envvar='DAMMIT_TEMP_DIR',
              hidden=True)
def show_directories_cmd(database_dir, temp_dir):
    '''List dammit directory locations.
    
    Locations come either from defaults or environment variables.'''

    print('Databases:', database_dir, file=sys.stderr)
    print('Temp root:', temp_dir, file=sys.stderr)
    print('Temp subdirs:', *(os.path.join(temp_dir, child) for child in TEMP_SUBDIRS),
          file=sys.stderr)


@config_group.command('clean-temp')
@click.option('--temp-dir',
              default=DEFAULT_TEMP_DIR,
              envvar='DAMMIT_TEMP_DIR')
@click.option('--envs', is_flag=True)
@click.option('--force', is_flag=True)
def clean_tmp_cmd(temp_dir, envs, force):
    ''' Clear out shared dammit temp files.
    '''

    temp_config_dir = os.path.join(temp_dir, WORKFLOW_CONFIG_TEMPDIR)
    print('Cleaning old config files in', temp_config_dir, file=sys.stderr)
    if not shutil.rmtree.avoids_symlink_attacks and not force:
        print('WARNING: platform and implementation is not restistant to symlink attacks.'
              ' Using the clean function could put this system at risk. Remove the directory'
              ' manually or use --force to preoceed.', file=sys.stderr)
    else:
        shutil.rmtree(temp_config_dir)
        if envs:
            conda_env_dir = os.path.join(temp_dir, CONDA_ENV_TEMPDIR)
            print('Removing conda environments in', conda_env_dir, file=sys.stderr)
            shutil.rmtree(conda_env_dir)

@config_group.command('show-default')
@click.pass_obj
@click.argument('config_file', type=click.Choice(['core', 'databases', 'pipelines']))
@click.option('--save', default='-', type=click.File('w'))
def show_default_cmd(defaults, config_file, save):
    ''' Show the selected default configuration file.'''

    save.write(yaml.dump(defaults.get(config_file)))


@config_group.command('busco-groups')
@click.pass_obj
def busco_groups_cmd(defaults):
    ''' Lists the available BUSCO group databases.'''

    db = defaults['databases']
    click.echo(' '.join(db['BUSCO']['groups'].keys()))
