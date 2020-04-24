#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2020
# File   : run.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 23.04.2020

import os
import sys

import click

from ..cli import component, CONFIG, ShortChoice


@component.group()
@click.pass_obj
@click.option('--database-dir',
              default=os.path.join(os.environ['HOME'], '.dammit', 'databases'),
              envvar='DAMMIT_DB_DIR',
              help='Directory to store databases. Existing'\
                    ' databases will not be overwritten.'\
                    ' By default, the database directory is'\
                    ' $HOME/.dammit/databases.')
@click.option('--busco-groups',
              default=['metazoa'],
              multiple=True,
              type=ShortChoice(list(CONFIG.databases['BUSCO']['groups'].keys()), case_sensitive=False),
              help='BUSCO group(s) to use/install.')
@click.option('--n_threads',
              type=int,
              help='Number of threads for overall workflow execution')
@click.option('--config-file',
              help='A YAML or JSON file providing values to override'\
                   ' built-in config. Advanced use only!')
@click.option('--busco-config-file',
              help='Path to an alternative BUSCO config'\
                   ' file; otherwise, BUSCO will attempt'\
                   ' to use its default installation'\
                   ' which will likely only work on'\
                   ' bioconda. Advanced use only!')
@click.option('--pipeline',
             type=click.Choice(["default", "quick", "full", "nr"], case_sensitive=False),
             help='Which pipeline to use. Pipeline options:'\
                  ' quick: excludes: '
                  ' the Infernal Rfam tasks, the HMMER'\
                  ' Pfam tasks, and the LAST OrthoDB'\
                  ' and uniref90 tasks. Best for users'\
                  ' just looking to get basic stats'\
                  ' and conditional reciprocal best'\
                  ' LAST from a protein database.'\
                  ' \nfull: '\
                  'Run a "complete" annotation; includes'\
                  ' uniref90, which is left out of the'\
                  ' default pipeline because it is huge'\
                  ' and homology searches take a long'\
                  ' time.'
                  ' nr: '\
                  ' Also include annotation to NR database, which'\
                  ' is left out of the default and "full"'\
                  ' pipelines because it is huge and'\
                  ' homology searches take a long time.'\
                  ' More info  at https://dib-lab.github.io/dammit.')
def run(config,
        database_dir,
        busco_groups, 
        n_threads, 
        config_file, 
        busco_config_file, 
        pipeline):

    if config_file:
        with open(config_file) as fp:
            config.core.update(yaml.safe_load(fp))
    
    if database_dir:
        config.core['database_dir'] = database_dir
    config.core['busco_groups'] = busco_groups

    if not n_threads and 'n_threads' not in config.core:
        config.core['n_threads'] = 1
    elif n_threads:
        config.core['n_threads'] = n_threads

    config.core['busco_config_file'] = busco_config_file

    if not pipeline and 'pipeline' not in config.core:
        config.core['pipeline'] = 'default'
    if pipeline:
        config.core['pipeline'] = pipeline

    click.echo(database_dir)


@run.command()
@click.pass_obj
def annotate(config):
    ''' The main annotation pipeline. Calculates assembly stats;
    runs BUSCO; runs LAST against OrthoDB (and optionally uniref90),
    HMMER against Pfam, Inferal against Rfam, and Conditional Reciprocal
    Best-hit Blast against user databases; and aggregates all results in
    a properly formatted GFF3 file.'''
    pass


@run.command()
@click.pass_obj
@click.option('--install', is_flag=True,
              help='Install missing databases. Downloads'
                   ' and preps where necessary')
def databases(config, install):
    pass
