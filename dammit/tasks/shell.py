# Copyright (C) 2015-2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the BSD license.  See the LICENSE file for details.

import os
import sys
import hashlib
import gzip


from doit.action import CmdAction
from doit.exceptions import TaskFailed
from doit.tools import LongRunning, run_once
from doit.task import clean_targets

from dammit.tasks.utils import clean_folder
from dammit.utils import which, doit_task


def hashfile(path, hasher=None, blocksize=65536):
    """
    A function to hash files.

    See: http://stackoverflow.com/questions/3431825
    """

    if hasher is None: hasher = hashlib.md5()

    try:
        try:
            f = gzip.open(path, "rb")
            buf = f.read(blocksize)
        except OSError:
            f = open(path, "rb")
            buf = f.read(blocksize)
    except FileNotFoundError:
        raise RuntimeError('Function hashfile could not find referenced file.'\
                           ' Is there a problem with curl?')

    while len(buf) > 0:
        hasher.update(buf)
        buf = f.read(blocksize)

    f.close()
    return hasher.hexdigest()


def check_hash(target_fn, expected):
    print('    * Checking hash of {0}'.format(target_fn), file=sys.stderr)
    if expected == hashfile(target_fn):
        return True
    else:
        os.remove(target_fn)
        return TaskFailed('{0} has non-matching hash; download error?'.format(target_fn))


@doit_task
def get_download_task(url, target_fn, md5=None, metalink=None):
    '''Creates a doit task to download the given URL.
    
    Args:
        url (str): URL to download.
        target_fn (str): Target for the download.
    Returns:
        dict: doit task.
    '''

    cmd = ['curl', '-o', target_fn]
    if metalink is not None:
        cmd.extend(['--metalink', metalink])
    cmd.append(url)
    cmd = ' '.join(cmd)
    name = 'download:{0}'.format(os.path.basename(target_fn))

    actions = [LongRunning(cmd)]

    if md5 is not None:
        actions.append((check_hash, [target_fn, md5]))

    return {'name': name,
            'actions': actions,
            'targets': [target_fn],
            'clean': [clean_targets],
            'uptodate': [True]}


@doit_task
def get_untargz_task(archive_fn, target_dir, label=None):
    '''Create a doit task to untar and gunip a *.tar.gz archive.

    Args:
        archive_fn (str): The .tar.gz file.
        target_dir (str): The folder to untar into.
        label (str): Optional label to resolve doit task name conflicts.
    Returns:
        dict: doit task.
    '''

    if label is None:
        label = os.path.basename(url)

    cmd = 'tar -xzf -C {target_dir} {archive_fn}'.format(**locals())
    name = 'untargz:{0}-{1}'.format(os.path.basename(target_dir), label)
    done = os.path.join(target_dir, name) + '.done'
    touch = 'touch {done}'.format(done=done)

    return {'name': name,
            'actions': ['mkdir -p {0}'.format(target_dir),
                        LongRunning(cmd),
                        touch],
            'targets': [done],
            'clean': [(clean_folder, [target_dir])],
            'uptodate': [True]}


@doit_task
def get_gunzip_task(archive_fn, target_fn):
    '''Create a doit task to gunzip a gzip archive.

    Args:
        archive_fn (str): The gzip file.
        target_fn (str): Output filename.
    Returns:
        dict: doit task.
    '''

    name = 'gunzip:{0}'.format(os.path.basename(target_fn))
    cmd = 'gunzip -c {archive_fn} > {target_fn}'.format(**locals())

    return {'name': name,
            'actions': [LongRunning(cmd)],
            'file_dep': [archive_fn],
            'targets': [target_fn],
            'clean': [clean_targets],
            'uptodate': [True]}


@doit_task
def get_download_and_gunzip_task(url, target_fn):
    '''Create a doit task which downloads and gunzips a file.

    Args:
        url (str): URL to download.
        target_fn (str): Target file for the download.
    Returns:
        dict: doit task.
    '''
    cmd = 'curl {url} | gunzip -c > {target_fn}'.format(**locals())

    name = 'download_and_gunzip:{0}'.format(os.path.basename(target_fn))

    return {'name': name,
            'actions': [LongRunning(cmd)],
            'targets': [target_fn],
            'clean': [clean_targets],
            'uptodate': [True]}


@doit_task
def get_download_and_untar_task(url, target_dir, label=None):
    '''Create a doit task to download a file and untar it in the
    given directory.

    Args:
        url (str): URL to download.
        target_dir (str: Directory to put the untarred folder in.
        label (str): Optional label to resolve doit name conflicts when putting
                     multiple results in the same folder.
    Returns:
        dict: doit task.
    '''

    if label is None:
        label = os.path.basename(url)

    cmd1 = 'mkdir -p {target_dir}; curl {url} | tar -xz -C {target_dir}'.format(**locals())
    name = 'download_and_untar:{0}-{1}'.format(os.path.basename(target_dir), label)
    done = os.path.join(target_dir, name) + '.done'
    cmd2 = 'touch {done}'.format(done=done)

    return {'name': name,
            'actions': [LongRunning(cmd1), cmd2],
            'targets': [done],
            'clean': [(clean_folder, [target_dir])],
            'uptodate': [True]}


@doit_task
def get_cat_task(file_list, target_fn):
    '''Create a doit task to `cat` together the given files and pipe the
    result to the given target.

    Args:
        file_list (list): The files to `cat`.
        target_fn (str): The target file.

    Returns:
        dict: A doit task.
    '''

    cmd = 'cat {files} > {t}'.format(files=' '.join(file_list), t=target_fn)

    return {'name': 'cat:' + os.path.basename(target_fn),
            'actions': [cmd],
            'file_dep': file_list,
            'targets': [target_fn],
            'clean': [clean_targets]}


@doit_task
def get_link_file_task(src, dst=''):
    ''' Soft-link file to the current directory, or to the destination
    target if given.

    Args:
        src (str): The file to link.
        dst (str): The destination; by default, the current directory.

    Returns:
        dict: A doit task.
    '''
    cmd = 'ln -fs {src} {dst}'.format(src=src, dst=dst)
    return {'name': 'ln:' + os.path.basename(src) + ('-' + dst if dst else ''),
            'actions': [cmd],
            'file_dep': [src],
            'targets': [os.path.basename(src) if not dst else dst],
            'uptodate': [run_once],
            'clean': [clean_targets]}
