#!/usr/bin/env python
import os

from doit.action import CmdAction
from doit.tools import title_with_actions, LongRunning
from doit.task import clean_targets

from .utils import which, doit_task, convert_pathlib
from .parallel import parallel_fasta, multinode_parallel_fasta, get_filesize_task


@doit_task
@convert_pathlib
def get_lastdb_task(db_fn, db_out_prefix=None, prot=True, params=None,
                    task_dep=None):
    '''Create a pydoit task to run lastdb.

    WARNING: This does not define a file_dep, to make sure it doesn't
    get executed when the dependency and targets already exist. This means
    that if the db_fn is created by another task, it MUST be defined before the
    lastdb task. This is a kludge which will eventually be fixed...

    Args:
        db_fn (str): The FASTA file to format.
        db_out_prefix (str): Prefix for the database files. Same as db_fn
                             if None (default).
        prot (bool): True if a protein FASTA, False otherwise.
        params (list): A list of additional parameters.
    Returns:
        dict: A pydoit task.
    '''

    cmd = [which('lastdb')]
    if prot:
        cmd.append('-p')
    if params is not None:
        cmd.extend([str(p) for p in params])
    if db_out_prefix is None:
        db_out_prefix = db_fn
    cmd.extend([db_out_prefix, db_fn])
    cmd = ' '.join(cmd)

    name = 'lastdb:' + os.path.basename(db_out_prefix)

    task_d =  {'name': name,
              'title': title_with_actions,
              'actions': [cmd],
              'targets': ['{0}.prj'.format(db_out_prefix)],
              'uptodate': [True],
              'clean': [clean_targets]}
    if task_dep is not None:
        task_d['task_dep'] = task_dep

    return task_d


@doit_task
@convert_pathlib
def get_lastal_task(query, db, out_fn, translate=False, frameshift=15,
                    cutoff=0.00001, n_threads=1, n_nodes=None, params=None):
    '''Create a pydoit task to run lastal

    Args:
        query (str): The file with the query sequences.
        db (str): The database file prefix.
        out_fn (str): Destination file for alignments.
        translate (bool): True if query is a nucleotide FASTA.
        frameshift (int): Frameshift penalty for translated alignment.
        n_threads (int): Number of threads to run with.
        n_nodes (int): Number of nodes for cluster environments. n_threads
                       will be the processors per node.
    Returns:
        dict: A pydoit task.
    '''

    lastal_exc = which('lastal')
    name = 'lastal:{0}'.format(os.path.join(out_fn))

    def lastal_cmd(file_size, cutoff=cutoff):
        cmd = [lastal_exc]
        if translate:
            cmd.append('-F' + str('frameshift'))
        if cutoff is not None:
            cutoff = round(1.0 / cutoff, 2)
            cmd.append('-D' + str(cutoff))
        if params is not None:
            cmd.extend([str(p) for p in params])
        cmd.append(db)
        cmd = '"{0}"'.format(' '.join(cmd))

        if n_nodes is None:
            parallel = parallel_fasta(query, n_threads, file_size)
        else:
            parallel = multinode_parallel_fasta(query, n_threads, n_nodes,
                                                file_size)

        cmd = [parallel, cmd, '<', query, '>', out_fn]
        return ' '.join(cmd)

    return {'name': name,
            'title': title_with_actions,
            'actions': [LongRunning(lastal_cmd)],
            'targets': [out_fn],
            'file_dep': [db + '.prj'],
            'getargs': {'file_size': ('get_filesize:{0}'.format(query),
                                      'size')},
            'clean': [clean_targets]}



def lastal(query, db, out_fn, **lastal_kwds):

    if cfg is None:
        cfg = config['lastal']
        
    yield get_filesize_task(query)
    yield get_lastal_task(query, db, out_fn, **lastal_kwds)
