#/usr/bin/env python
import os

from doit.action import CmdAction
from doit.tools import title_with_actions
from doit.task import clean_targets

from .common import which
from .tasks import doit_task
from .parallel import parallel_fasta, multinode_parallel_fasta


@doit_task
def get_lastdb_task(db_fn, db_out_prefix, lastdb_cfg, prot=True):
    '''Create a pydoit task to run lastdb.

    WARNING: This does not define a file_dep, to make sure it doesn't
    get executed when the dependency and targets already exist. This means
    that if a task acquires the database, it MUST be defined before the
    lastdb task.

    Args:
        db_fn (str): The FASTA file to format.
        db_out_prefix (str): Prefix for the database files.
        lastdb_cfg (dict): Config for the command. Shoud contain an entry
        named "params" storing a str.
        prot (bool): True if a protein FASTA, False otherwise.
    Returns:
        dict: A pydoit task.
    '''

    exc = which('lastdb')
    params = lastdb_cfg['params']
    if prot:
        params += ' -p'

    cmd = '{exc} {params} {db_out_prefix} {db_fn}'.format(**locals())

    name = 'lastdb:' + os.path.basename(db_out_prefix)

    return {'name': name,
            'title': title_with_actions,
            'actions': [cmd],
            'targets': ['{0}.prj'.format(db_out_prefix)],
            'uptodate': [True],
            'clean': [clean_targets]}


@doit_task
def get_lastal_task(query, db, out_fn, cfg, translate=False, 
                    cutoff=0.00001, n_threads=1, n_nodes=None):
    '''Create a pydoit task to run lastal

    Args:
        query (str): The file with the query sequences.
        db (str): The database file prefix.
        out_fn (str): Destination file for alignments.
        translate (bool): True if query is a nucleotide FASTA.
        n_threads (int): Number of threads to run with.
        cfg (dict): Config, must contain key params holding str.
    Returns:
        dict: A pydoit task.
    '''

    lastal_exc = which('lastal')

    params = cfg['params']
    lastal_cmd = [lastal_exc]
    if translate:
        lastal_cmd.append('-F' + str(cfg['frameshift']))
    if cutoff is not None:
        cutoff = round(1.0 / cutoff, 2)
        lastal_cmd.append('-D' + str(cutoff))
    lastal_cmd.append(db)
    lastal_cmd = '"{0}"'.format(' '.join(lastal_cmd))

    if n_nodes is None:
        parallel = parallel_fasta(query, n_threads)
    else:
        parallel = multinode_parallel_fasta(query, n_threads, n_nodes)

    cmd = [parallel, lastal_cmd, '<', query, '>', out_fn]
    cmd = ' '.join(cmd)

    name = 'lastal:{0}'.format(os.path.join(out_fn))

    return {'name': name,
            'title': title_with_actions,
            'actions': [cmd],
            'targets': [out_fn],
            'file_dep': [db + '.prj'],
            'clean': [clean_targets]}



