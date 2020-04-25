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

from ..config import CONFIG
from ..meta import __path__
from ..utils import ShortChoice



@click.group('run')
@click.pass_obj
@click.option('--database-dir',
              default=os.path.join(os.environ['HOME'], '.dammit', 'databases'),
              envvar='DAMMIT_DB_DIR',
              help='Directory to store databases. Existing'\
                    ' databases will not be overwritten.'\
                    ' By default, the database directory is'\
                    ' $HOME/.dammit/databases.')
@click.option('--busco-group',
              default=['metazoa'],
              multiple=True,
              type=ShortChoice(list(CONFIG.databases['BUSCO']['groups'].keys()), case_sensitive=False),
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
              busco_group, 
              n_threads, 
              config_file, 
              busco_config_file, 
              pipeline):
    ''' Run the annotation pipeline or install databases.
    '''

    if config_file:
        with open(config_file) as fp:
            config.core.update(yaml.safe_load(fp))
    
    if database_dir:
        config.core['database_dir'] = database_dir
    config.core['busco_groups'] = busco_group

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
def annotate_cmd(config,
                 transcriptome,
                 name,
                 evalue,
                 output_dir,
                 user_database):
    ''' The main annotation pipeline. Calculates assembly stats;
    runs BUSCO; runs LAST against OrthoDB (and optionally uniref90),
    HMMER against Pfam, Inferal against Rfam, and Conditional Reciprocal
    Best-hit Blast against user databases; and aggregates all results in
    a properly formatted GFF3 file.'''
    
    print(transcriptome, name, evalue, output_dir, user_database)


@run_group.command('databases')
@click.pass_obj
@click.option('--install', is_flag=True,
              help='Install missing databases. Downloads'
                   ' and preps where necessary')
def databases_cmd(config, install):
    ''' The database preparation pipeline.
    '''
    workflow_file = os.path.join(__path__, 'workflows', 'dammit.snakefile')
    config_file = os.path.join(__path__, 'config.yml')
    database_dir = config.core['database_dir']
    pipeline_config = config.pipelines['pipelines'][config.core['pipeline']]
    n_threads = config.core['n_threads']

    cmd = ["snakemake", "-s", workflow_file, "--configfile", config_file]

    targets = generate_database_targets(pipeline_config, database_dir, config)
    
    config = ["--config", f"db_dir={database_dir}"]
    cmd.extend(config)

    helpful_args = ["-p", "--nolock", "--use-conda", "--rerun-incomplete", "-k", "--cores", f"{n_threads}"]
    cmd.extend(helpful_args)

    # better way to do this?
    if not install:
        cmd.append("--dry-run")

    # finally, add targets
    cmd.extend(targets)
    print("Command: " + " ".join(cmd))

    try:
        subprocess.check_call(cmd)
    except subprocess.CalledProcessError as e:
        print(f'Error in snakemake invocation: {e}', file=sys.stderr)
        return e.returncode


def generate_database_targets(pipeline_info, database_dir, config):
    targets = []
    pipeline_databases = pipeline_info["databases"]

    # TODO: fix BUSCO
    if "BUSCO" in pipeline_databases:
        pipeline_databases.remove("BUSCO")

    for db in pipeline_databases:
        fn = config.databases[db]["filename"]
        #out_suffixes = [""] # testing: ONLY DOWNLOAD
        out_suffixes = config.databases[db]["output_suffix"]
        targets += [fn + suffix for suffix in out_suffixes]
    targets = [os.path.join(database_dir, targ) for targ in targets]

    return targets

'''
def generate_annotation_targets
    # generate annotation targets
    if annot:
        #
        # Generating targets for annotate subcommand
        #
        transcriptome = config.core['transcriptome']
        if any([transcriptome.endswith(".fa"),
                transcriptome.endswith(".fasta")]):
            transcriptome_name = os.path.basename(transcriptome).rsplit(".fa")[0]
        else:
            raise ValueError('input transcriptome file must end with ".fa" or ".fasta"')

        if self.args.output_dir:
            out_dir = self.args.output_dir
        else:
            out_dir = transcriptome_name + self.config_d["dammit_dir"]

        annotation_programs = pipeline_info["programs"]
        annotation_databases = pipeline_info["databases"]
        output_suffixes = []
        # not complete yet. need to include database name in annotation targ, where relevant
        # not sure how to represent this in the config.yml. databases arg for prog?
        for prog in annotation_programs:
            prog_suffixes = self.config_d[prog]["output_suffix"]
            prog_databases = self.config_d[prog].get("databases")
            if prog_databases:
                # only consider databases we're running in this pipeline
                dbs_to_add = [db for db in prog_databases if db in annotation_databases]
                # expand __database__ with appropriate databases
                db_suffixes = []
                for suffix in prog_suffixes:
                    if "__database__" in suffix:
                        for db in dbs_to_add:
                            db_suffixes.append(suffix.replace("__database__", db))
                    else:
                        db_suffixes.append(suffix)
                prog_suffixes = db_suffixes
            output_suffixes.extend(prog_suffixes)
        annotate_targets = [os.path.join(out_dir, transcriptome_name + suffix) for suffix in output_suffixes]
        targets+=annotate_targets
    return targets, db_dir, out_dir
'''
