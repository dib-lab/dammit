#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2020
# File   : config.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 23.04.2020

import os
import sys
import yaml

import click

from ..config import (parse_config, DEFAULT_DATABASES_DIR, DEFAULT_TEMP_DIR, TEMP_SUBDIRS)


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
def clean_tmp_cmd(temp_dir):
    ''' Clear out shared dammit temp files.
    '''
    print('Cleaning files from:', *(os.path.join(temp_dir, child) for child in TEMP_SUBDIRS),
          file=sys.stderr)


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
