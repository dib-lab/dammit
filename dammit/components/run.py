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

import click
import psutil

from rich import print as richprint
from rich.console import RenderGroup
from rich.padding import Padding
from rich.panel import Panel as RichPanel
from rich.tree import Tree as RichTree

from ..config import DATABASES
from ..meta import __path__, __time__
from ..utils import ShortChoice, write_yaml, update_nested_dict, create_dirs


@click.group('run')
@click.pass_obj
@click.option('--database-dir',
              envvar='DAMMIT_DB_DIR',
              help='Directory to store databases. Existing'\
                    ' databases will not be overwritten.')
@click.option('--conda-dir',
              envvar='DAMMIT_CONDA_DIR',
              help='Directory to store snakemake-created conda environments.')
@click.option('--temp-dir',
              envvar='DAMMIT_TEMP_DIR',
              help='Directory to store dammit temp files.')
@click.option('--busco-group',
              multiple=True,
              type=ShortChoice(list(DATABASES['busco']['lineages']), case_sensitive=False),
              help='BUSCO group(s) to use/install.')
@click.option('--n-threads',
              envvar='DAMMIT_N_THREADS',
              type=int,
              help='Number of threads for overall workflow execution')
@click.option('--max-threads-per-task',
               type=int,
               help='Max threads to use for a single step.')
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
              conda_dir,
              temp_dir,
              busco_group,
              n_threads,
              max_threads_per_task,
              busco_config_file,
              pipeline):
    ''' Run the annotation pipeline or install databases.
    '''

    print(config.banner, file=sys.stderr)
    config.gui.param_tree = RichTree('📔 Configuration')


    if database_dir:
        config.core['database_dir'] = os.path.abspath(database_dir)

    if temp_dir:
        config.core['temp_dir'] = os.path.abspath(temp_dir)
    
    if conda_dir:
        config.core['conda_env_dir'] = os.path.abspath(conda_dir)

    create_dirs([config.core[k] for k in ['database_dir', 'temp_dir', 'conda_env_dir']])

    if busco_group:
        config.core['busco_groups'] = list(busco_group)

    if not n_threads or n_threads == 0:
        config.core['n_threads'] = psutil.cpu_count(logical=False)
    else:
        config.core['n_threads'] = n_threads
    
    if not max_threads_per_task:
        config.core['max_threads_per_task'] = config.core['n_threads']
    else:
        config.core['max_threads_per_task'] = min(config.core['n_threads'], max_threads_per_task)

    if busco_config_file:
        config.core['busco']['configfile'] = busco_config_file
    else:
        config.core['busco']['configfile'] = os.path.join(__path__, config.core["busco"]["configfile"])

    if pipeline:
        config.core['pipeline'] = pipeline

    config.gui.core_params = config.gui.param_tree.add(
        RenderGroup(
            '🌐 Global params',
            RichPanel(f"{'Pipeline:'.ljust(25)} {config.core['pipeline']}\n"
                      f"{'BUSCO groups:'.ljust(25)} {', '.join(config.core['busco_groups'])}\n"
                      f"{'Databases:'.ljust(25)} {config.core['database_dir']}\n"
                      f"{'Conda Environments:'.ljust(25)} {config.core['conda_env_dir']}\n"
                      f"{'Threads (total):'.ljust(25)} {config.core['n_threads']}\n"
                      f"{'Threads (per-task):'.ljust(25)} {config.core['max_threads_per_task']}",
                      expand=True),
        )
    )


@run_group.command('annotate', context_settings=dict(
    ignore_unknown_options=True
))
@click.pass_obj
@click.argument('transcriptome')
@click.option('-n', '--base-name',
              help='Base name to use for renaming the'\
                   ' input transcripts. The new names'\
                   ' will be of the form <name>_<X>.'\
                   ' It should not have spaces, pipes,'\
                   ' ampersands, or other characters'\
                   ' with special meaning to BASH.'\
                   ' Superseded by --regex-rename.')
@click.option('--regex-rename',
              help='Rename transcripts using a regex pattern. The regex should follow '
                   ' Python `re` format and contain a named field keyed'\
                   ' as `name` that extracts the desired string. For example, providing'\
                   ' (?P<name>^[a-zA-Z0-9]+) will match from the beginning of the sequence header'\
                   ' up to the first non-alphanumeric symbol.'\
                   ' Supersedes --base-name.')
@click.option('--rename/--no-rename', default=None,
               help='If --no-rename, original transcript names are preserved'\
                    ' in the final annotated FASTA. --base-name is'\
                    ' still used in intermediate files. If --rename (the default '\
                    ' behavior), the renamed transcript names are used in the final '\
                    ' annotated FASTA.')
@click.option('-e', '--global-evalue',
              type=float,
              help='global e-value cutoff for similarity searches.')
@click.option('-o', '--output-dir',
              help='Output directory. By default this will'\
                   ' be the name of the transcriptome file'\
                   ' with `.dammit` appended')
@click.option('-u', '--user-database',
              multiple=True,
              help='Optional additional protein databases. '\
                   ' These will be searched with CRB-blast.')
@click.option('--dry-run', is_flag=True)
@click.argument('extra_snakemake_args', nargs=-1, type=click.UNPROCESSED)
def annotate_cmd(config,
                 transcriptome,
                 base_name,
                 regex_rename,
                 rename,
                 global_evalue,
                 output_dir,
                 user_database,
                 dry_run,
                 extra_snakemake_args):
    ''' The main annotation pipeline. Calculates assembly stats;
    runs BUSCO; runs LAST against OrthoDB (and optionally uniref90),
    HMMER against Pfam, Inferal against Rfam, and Conditional Reciprocal
    Best-hit Blast against user databases; and aggregates all results in
    a properly formatted GFF3 file.'''

    config.core['command'] = 'annotate'

    # Strip the extension from the txome and use as the global prefix name
    if any([transcriptome.endswith(".fa"),
            transcriptome.endswith(".fasta")]):
        transcriptome_name = os.path.basename(transcriptome).rsplit(".fa")[0]
    else:
        raise ValueError('input transcriptome file must end with ".fa" or ".fasta"')
    config.core['transcriptome_name'] = transcriptome_name

    # if an output folder is not provided, use the default suffix
    if not output_dir:
        output_dir = os.path.abspath(transcriptome_name + config.core["dammit_dir_suffix"])
    else:
        output_dir = os.path.abspath(output_dir)
    config.core['output_dir'] = output_dir

    # we want the absolute path to the input txome
    config.core['input_transcriptome'] = os.path.abspath(transcriptome)

    if base_name:
        config.core['basename'] = base_name
    
    if regex_rename:
        config.core['regex_rename'] = regex_rename
    
    if rename is not None:
        config.core['rename'] = rename
    
    # the default global_evalue is null 
    if global_evalue:
        config.core['global_evalue'] = global_evalue

    config.core['user_dbs'] = wrangle_user_databases(user_database)

    # extract the given pipeline from the config
    pipeline_config = config.pipelines['pipelines'][config.core['pipeline']]

    # Wrangle snakemake files
    workflow_file = os.path.join(__path__, 'workflows', 'dammit.snakefile')

    # Create the output directory
    os.makedirs(output_dir, exist_ok=True)

    targets = generate_annotation_targets(pipeline_config, config)

    # config file for *this run*. gets put in the dammit output directory
    workflow_config_file = os.path.join(output_dir, 'run.config.yml')
    write_yaml(config.core, workflow_config_file)

    # Build snakemake command
    cmd = ["snakemake", "-s", workflow_file,
           "--configfiles", workflow_config_file,
           "--directory", output_dir]
    cmd.extend(snakemake_common_args(config.core['n_threads'],
                                     config.core['conda_env_dir']))
    if dry_run:
        cmd.append('--dry-run')

    cmd.extend(list(extra_snakemake_args))
    cmd.extend(targets)

    config.gui.annot_params = config.gui.param_tree.add(
        RenderGroup(
            '🖊️ Annotation params',
            RichPanel(f"{'Input:'.ljust(25)} {config.core['input_transcriptome']}\n"
                      f"{'Output:'.ljust(25)} {config.core['output_dir']}\n"
                      f"{'E-value Cutoff (global):'.ljust(25)} {config.core['global_evalue']}\n"
                      f"{'User databases:'.ljust(25)} {', '.join(config.core['user_dbs'])}\n"
                      f"{'Run config:'.ljust(25)} {workflow_config_file}\n"
                      f"{'Extra Snakemake args:'.ljust(25)} {extra_snakemake_args}",
                      expand=True),
        )
    )

    config.gui.snakemake_cmd = config.gui.param_tree.add(
        RenderGroup(
            '🐍 Snakemake Invocation',
            RichPanel(" ".join(cmd))
        )
    )
    richprint(config.gui.param_tree, file=sys.stderr)

    richprint('\n▶️  Beginning workflow execution...\n', file=sys.stderr)

    try:
        subprocess.check_call(cmd)
    except subprocess.CalledProcessError as e:
        print(f'Error in snakemake invocation!', file=sys.stderr)
        print('If the error is a MissingInputException, you probably need to install '
              'the dammit databases. Run `dammit run databases --help` for more information.', file=sys.stderr)
        sys.exit(e.returncode)


@run_group.command('databases', context_settings=dict(
    ignore_unknown_options=True
))
@click.pass_obj
@click.option('--install', is_flag=True,
              help='Install missing databases. Downloads'
                   ' and preps where necessary.')
@click.argument('extra_snakemake_args', nargs=-1, type=click.UNPROCESSED)
def databases_cmd(config, install, extra_snakemake_args):
    ''' The database preparation pipeline. The database installation directory is set
    after the `run` subcommand (invoke `dammit run --help` for more information). You
    can also set the ${DAMMIT_DB_DIR} environment variable. For example:

        dammit run --database-dir /path/to/database/dir databases --install
    '''
    config.core['command'] = 'databases'

    workflow_file = os.path.join(__path__, 'workflows', 'dammit.snakefile')
    
    output_dir = os.path.join(config.core['temp_dir'], 
                               f'run.databases.{__time__}')
    config.core['output_dir'] = output_dir
    os.makedirs(output_dir, exist_ok=True)

    workflow_config_file = os.path.join(output_dir,
                                        f'config.yml')
    write_yaml(config.core, workflow_config_file)

    pipeline_config = config.pipelines['pipelines'][config.core['pipeline']]

    cmd = ["snakemake", "-s", workflow_file, "--configfiles", workflow_config_file]

    targets = generate_database_targets(pipeline_config, config)
    cmd.extend(snakemake_common_args(config.core['n_threads'],
                                     config.core['conda_env_dir']))
    cmd.extend(list(extra_snakemake_args))

    # better way to do this?
    if not install:
        cmd.append("--dry-run")
        cmd.append("--detailed-summary")
        #cmd.append("--list-target-rules")

    # finally, add targets
    cmd.extend(targets)

    config.gui.annot_params = config.gui.param_tree.add(
        RenderGroup(
            '🗄️ Database params',
            RichPanel(f"{'Run config:'.ljust(25)} {workflow_config_file}\n"
                      f"{'Extra Snakemake args:'.ljust(25)} {extra_snakemake_args}",
                      expand=True),
        )
    )

    config.gui.snakemake_cmd = config.gui.param_tree.add(
        RenderGroup(
            '🐍 Snakemake Invocation',
            RichPanel(" ".join(cmd))
        )
    )
    richprint(config.gui.param_tree, file=sys.stderr)

    if install:
        richprint('\n▶️  Beginning workflow execution...\n', file=sys.stderr)
    else:
        richprint('\n❔ Database status (use with `--install` to run database installation pipeline):\n',
                  file=sys.stderr)

    try:
        subprocess.check_call(cmd)
    except subprocess.CalledProcessError as e:
        print(f'Error in snakemake invocation: {e}', file=sys.stderr)
        sys.exit(e.returncode)
    

def generate_database_targets(pipeline_info, config):
    targets = []
    pipeline_databases = pipeline_info["databases"]
    database_dir = config.core['database_dir']

    for db in pipeline_databases:
        try:
            fn = config.databases[db]["filename"]
        except:
            continue
        out_suffixes = config.databases[db]["output_suffix"]
        targets += [fn + suffix for suffix in out_suffixes]
    targets = [os.path.join(database_dir, targ) for targ in targets]

    return targets


def generate_annotation_targets(pipeline_info, config):
    annotation_programs = pipeline_info["programs"]
    annotation_databases = pipeline_info.get("databases", [])
    transcriptome_name = config.core['transcriptome_name']
    output_dir = config.core['output_dir']
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


def snakemake_common_args(n_threads, conda_env_dir):
    args = ["-p", "--nolock", "--conda-frontend", "mamba",
            "--use-conda", "--conda-prefix", conda_env_dir,
            "--rerun-incomplete", "-k", "--cores", str(n_threads)]
    return args


def wrangle_user_databases(user_dbs):
    dbs = {}
    for udb in user_dbs:
        udb_path = os.path.abspath(os.path.expanduser(udb)) # get absolute path, expanding any `~`
        udb_name = os.path.basename(udb_path)
        dbs[udb_name] = udb_path
    return dbs
