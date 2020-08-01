#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2020
# File   : run.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 23.04.2020

import os
import subprocess
import sys
import yaml

import click

from ..config import (CONFIG, DEFAULT_CONFIG_DIR, CONDA_ENV_TEMPDIR,
                      WORKFLOW_CONFIG_TEMPDIR, DEFAULT_DATABASES_DIR,
                      DEFAULT_TEMP_DIR, create_tempdirs)
from ..meta import __path__, __time__
from ..utils import ShortChoice, read_yaml, write_yaml, update_nested_dict



@click.group('run')
@click.pass_obj
@click.option('--database-dir',
              default=DEFAULT_DATABASES_DIR,
              envvar='DAMMIT_DB_DIR',
              help='Directory to store databases. Existing'\
                    ' databases will not be overwritten.')
@click.option('--temp-dir',
              default=DEFAULT_TEMP_DIR,
              envvar='DAMMIT_TEMP_DIR',
              help='Directory to store dammit temporary files.'\
                   ' These will include the final workflow configs'\
                   ' for individual database runs as well as conda'\
                   ' environments for snakemake rules.')
@click.option('--busco-group',
              default=['metazoa_odb10'],
              multiple=True,
              type=ShortChoice(list(CONFIG.databases['busco']['lineages']), case_sensitive=False),
              help='BUSCO group(s) to use/install.')
@click.option('--n-threads',
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
def run_group(config,
              database_dir,
              temp_dir,
              busco_group,
              n_threads,
              config_file,
              busco_config_file,
              pipeline):
    ''' Run the annotation pipeline or install databases.
    '''

    if config_file:
        user_config = read_yaml(config_file)
        update_nested_dict(config.core, user_config)

    config.core['db_dir'] = database_dir
    os.makedirs(database_dir, exist_ok=True)

    config.core['temp_dir'] = temp_dir
    create_tempdirs(temp_dir)

    config.core['busco_groups'] = list(busco_group)

    if not n_threads and 'n_threads' not in config.core:
        config.core['n_threads'] = 1
    elif n_threads:
        config.core['n_threads'] = n_threads

    if busco_config_file:
        config.core['busco_config_file'] = busco_config_file
    else:
        config.core['busco_config_file'] = os.path.join(database_dir, config.databases["busco"]["filename"])

    if not pipeline and 'pipeline' not in config.core:
        config.core['pipeline'] = 'default'
    if pipeline:
        config.core['pipeline'] = pipeline

    click.echo(database_dir)


@run_group.command('annotate')
@click.pass_obj
@click.argument('transcriptome')
@click.option('-n', '--name', default='Transcript',
              help='Base name to use for renaming the'\
                   ' input transcripts. The new names'\
                   ' will be of the form <name>_<X>.'\
                   ' It should not have spaces, pipes,'\
                   ' ampersands, or other characters'\
                   ' with special meaning to BASH.')
@click.option('-e', '--evalue',
              default=1e-5,
              type=float,
              help='e-value cutoff for similarity searches.')
@click.option('-o', '--output-dir',
              default=None,
              help='Output directory. By default this will'\
                   ' be the name of the transcriptome file'\
                   ' with `.dammit` appended')
@click.option('-u', '--user-database',
              multiple=True,
              help='Optional additional protein databases. '\
                   ' These will be searched with CRB-blast.')
@click.option('--dry-run', is_flag=True)
def annotate_cmd(config,
                 transcriptome,
                 name,
                 evalue,
                 output_dir,
                 user_database,
                 dry_run):
    ''' The main annotation pipeline. Calculates assembly stats;
    runs BUSCO; runs LAST against OrthoDB (and optionally uniref90),
    HMMER against Pfam, Inferal against Rfam, and Conditional Reciprocal
    Best-hit Blast against user databases; and aggregates all results in
    a properly formatted GFF3 file.'''

    # Handle config updates
    if any([transcriptome.endswith(".fa"),
            transcriptome.endswith(".fasta")]):
        transcriptome_name = os.path.basename(transcriptome).rsplit(".fa")[0]
    else:
        raise ValueError('input transcriptome file must end with ".fa" or ".fasta"')
    config.core['transcriptome_name'] = transcriptome_name

    if output_dir is None:
        output_dir = os.path.abspath(transcriptome_name + config.core["dammit_dir_suffix"])
    config.core['dammit_dir'] = output_dir
    config.core['input_transcriptome'] = os.path.abspath(transcriptome)

    config.core['user_dbs'] = wrangle_user_databases(user_database)

    pipeline_config = config.pipelines['pipelines'][config.core['pipeline']]

    # Wrangle snakemake files
    workflow_file = os.path.join(__path__, 'workflows', 'dammit.snakefile')

    # Create the output directory
    os.makedirs(output_dir, exist_ok=True)

    targets = generate_annotation_targets(pipeline_config, config)

    workflow_config_file = os.path.join(output_dir, 'run.config.yml')
    print(f'Writing full run config to {workflow_config_file}', file=sys.stderr)
    write_yaml(config.core, workflow_config_file)

    # Build snakemake command
    cmd = ["snakemake", "-s", workflow_file,
           "--configfiles", workflow_config_file,
           "--directory", output_dir]
    cmd.extend(snakemake_common_args(config.core['n_threads'], config.core['temp_dir']))
    if dry_run:
        cmd.append('--dry-run')

    cmd.extend(targets)

    print("Command: " + " ".join(cmd), file=sys.stderr)
    try:
        subprocess.check_call(cmd)
    except subprocess.CalledProcessError as e:
        print(f'Error in snakemake invocation: {e}', file=sys.stderr)
        return e.returncode


@run_group.command('databases')
@click.pass_obj
@click.option('--install', is_flag=True,
              help='Install missing databases. Downloads'
                   ' and preps where necessary')
def databases_cmd(config, install):
    ''' The database preparation pipeline.
    '''
    workflow_file = os.path.join(__path__, 'workflows', 'dammit.snakefile')
    database_dir = config.core['db_dir']

    workflow_config_file = os.path.join(config.core['temp_dir'],
                                        WORKFLOW_CONFIG_TEMPDIR,
                                        f'databases.config.{__time__}.yml')
    print(f'Writing full run config to {workflow_config_file}', file=sys.stderr)
    write_yaml(config.core, workflow_config_file)

    pipeline_config = config.pipelines['pipelines'][config.core['pipeline']]

    cmd = ["snakemake", "-s", workflow_file, "--configfiles", workflow_config_file]

    targets = generate_database_targets(pipeline_config, config)
    cmd.extend(snakemake_common_args(config.core['n_threads'], config.core['temp_dir']))

    # better way to do this?
    if not install:
        cmd.append("--dry-run")

    # finally, add targets
    cmd.extend(targets)
    print("Command: " + " ".join(cmd), file=sys.stderr)

    try:
        subprocess.check_call(cmd)
    except subprocess.CalledProcessError as e:
        print(f'Error in snakemake invocation: {e}', file=sys.stderr)
        return e.returncode


def generate_database_targets(pipeline_info, config):
    targets = []
    pipeline_databases = pipeline_info["databases"]
    database_dir = config.core['db_dir']

    for db in pipeline_databases:
        fn = config.databases[db]["filename"]
        out_suffixes = config.databases[db]["output_suffix"]
        targets += [fn + suffix for suffix in out_suffixes]
    targets = [os.path.join(database_dir, targ) for targ in targets]

    return targets


def generate_annotation_targets(pipeline_info, config):
    annotation_programs = pipeline_info["programs"]
    annotation_databases = pipeline_info.get("databases", [])
    transcriptome_name = config.core['transcriptome_name']
    output_dir = config.core['dammit_dir']
    user_dbs = config.core['user_dbs']
    busco_lineages = config.core['busco_groups']

    output_suffixes = []

    for prog in annotation_programs:
        prog_suffixes = config.core[prog]["output_suffix"]
        prog_databases = config.core[prog].get("databases", [])
        if prog_databases:
            # now handle other databases. Only consider databases we're running in this pipeline
            dbs_to_add = [db for db in prog_databases if db in annotation_databases]
            # expand __database__ with appropriate databases
            if prog_suffixes and '__userdatabase__' in prog_suffixes[0]:
                prog_suffixes = [suffix.replace("__userdatabase__", db) for db in user_dbs.keys() for suffix in prog_suffixes]
            if dbs_to_add:
                prog_suffixes = [suffix.replace("__database__", db) for db in dbs_to_add for suffix in prog_suffixes]
            if busco_lineages:
                prog_suffixes = [suffix.replace("__buscolineage__", db) for db in busco_lineages for suffix in prog_suffixes]
        output_suffixes.extend(prog_suffixes)
    targets = [os.path.join(output_dir, transcriptome_name + suffix) for suffix in output_suffixes]

    gff_files = [x for x in targets if x.endswith(".gff3")]
    config.core["gff_files"] = gff_files
    targets+=[os.path.join(output_dir, transcriptome_name + suffix) for suffix in config.core["output_suffix"]]

    return targets


def snakemake_common_args(n_threads, temp_dir):
    args = ["-p", "--nolock", "--conda-frontend", "mamba",
            "--use-conda", "--conda-prefix", os.path.join(temp_dir, CONDA_ENV_TEMPDIR),
            "--rerun-incomplete", "-k", "--cores", str(n_threads)]
    return args


def wrangle_user_databases(user_dbs):
    dbs = {}
    for udb in user_dbs:
        udb_path = os.path.abspath(os.path.expanduser(udb)) # get absolute path, expanding any `~`
        udb_name = os.path.basename(udb_path)
        dbs[udb_name] = udb_path
    return dbs
