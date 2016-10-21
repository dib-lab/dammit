#!/usr/bin/env python
import os

from doit.action import CmdAction
from doit.tools import title_with_actions, LongRunning
from doit.task import clean_targets

from .utils import clean_folder
from ..profile import profile_task
from ..utils import which, doit_task


@doit_task
@profile_task
def get_transdecoder_orf_task(input_filename,params=None):

    name = 'TransDecoder.LongOrfs:' + os.path.basename(input_filename)

    exc = which('TransDecoder.LongOrfs')
    cmd = [exc, '-t', input_filename]
    if params is not None:
        cmd.extend(params)
    cmd = ' '.join(cmd)

    return {'name': name,
            'title': title_with_actions,
            'actions': [cmd],
            'file_dep': [input_filename],
            'targets': [input_filename + '.transdecoder_dir/longest_orfs.pep'],
            'clean': [(clean_folder, [input_filename + '.transdecoder_dir'])]}

@doit_task
@profile_task
def get_transdecoder_predict_task(input_filename, pfam_filename, params=None):

    name = 'TransDecoder.Predict:' + os.path.basename(input_filename)

    exc = which('TransDecoder.Predict')
    cmd = [exc, '-t', input_filename, '--retain_pfam_hits', pfam_filename]
    if params is not None:
        cmd.extend(params)
    cmd = ' '.join(cmd)

    return {'name': name,
            'title': title_with_actions,
            'actions': [cmd],
            'file_dep': [input_filename,
                         input_filename + '.transdecoder_dir/longest_orfs.pep',
                         pfam_filename],
            'targets': [input_filename + '.transdecoder' + ext \
                        for ext in ['.bed', '.cds', '.pep', '.gff3', '.mRNA']],
            'clean': [clean_targets,
                     (clean_folder, [input_filename + '.transdecoder_dir'])]}

