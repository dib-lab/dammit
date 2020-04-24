#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2020
# File   : config.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 23.04.2020

import yaml

import click

from ..cli import component, CONFIG
from ..config import parse_config


@component.group(name='config')
@click.pass_context
def config_group(ctx):
    ''' Show dammit configuration information.'''

    core, databases, pipelines = parse_config()
    ctx.obj = {'core': core, 'databases': databases, 'pipelines': pipelines}


@config_group.command()
@click.pass_obj
@click.argument('config_file', type=click.Choice(['core', 'databases', 'pipelines']))
@click.option('--save', default='-', type=click.File('w'))
def show_default(defaults, config_file, save):
    ''' Show the selected default configuration file.'''

    save.write(yaml.dump(defaults.get(config_file)))


@config_group.command()
@click.pass_obj
def busco_groups(defaults):
    ''' Lists the available BUSCO group databases.'''

    db = defaults['databases']
    click.echo(' '.join(db['BUSCO']['groups'].keys()))
