#!/usr/bin/env python
from __future__ import print_function
import os

from doit.action import CmdAction
from doit.tools import title_with_actions, LongRunning
from doit.task import clean_targets

from .utils import clean_folder
from ..utils import which, doit_task


@doit_task
def get_download_task(url, target_fn):
    '''Creates a doit task to download the given URL.
    
    Args:
        url (str): URL to download.
        target_fn (str): Target for the download.
    Returns:
        dict: doit task.
    '''

    cmd = 'curl -o {target_fn} {url}'.format(**locals())
    name = 'download_gunzip:{0}'.format(os.path.basename(target_fn))

    return {'title': title_with_actions,
            'name': name,
            'actions': [LongRunning(cmd)],
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

    return {'title': title_with_actions,
            'name': name,
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
            'title': title_with_actions,
            'actions': [LongRunning(cmd1), cmd2],
            'targets': [done],
            'clean': [(clean_folder, [target_dir])],
            'uptodate': [True]}


@doit_task
def get_create_folder_task(folder):

    name = 'create_folder:{folder}'.format(**locals())

    return {'title': title_with_actions,
            'name': name,
            'actions': [(create_folder, [folder])],
            'targets': [folder],
            'uptodate': [run_once],
            'clean': [clean_targets] }


@doit_task
def get_cat_task(file_list, target_fn):

    cmd = 'cat {files} > {t}'.format(files=' '.join(file_list), t=target_fn)

    return {'title': title_with_actions,
            'name': 'cat:' + os.path.basename(target_fn),
            'actions': [cmd],
            'file_dep': file_list,
            'targets': [target_fn],
            'clean': [clean_targets]}


@doit_task
def get_link_file_task(src, dst=''):
    ''' Soft-link file to the current directory
    '''
    cmd = 'ln -fs {src} {dst}'.format(src=src, dst=dst)
    return {'title': title_with_actions,
            'name': 'ln:' + os.path.basename(src) + ('-' + dst if dst else ''),
            'actions': [cmd],
            'file_dep': [src],
            'targets': [os.path.basename(src) if not dst else dst],
            'uptodate': [run_once],
            'clean': [clean_targets]}
