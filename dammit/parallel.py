#!/usr/bin/env python
from __future__ import print_function

from .utils import which, doit_task

import os


def parallel_fasta(input_filename, output_filename, command, n_jobs, pbs=False):

    exc = which('parallel')
    cmd = ['cat', input_filename, '|', exc, '--progress', '--pipe', '-L', 2, '-N', 400,
           '--gnu', '-j', n_jobs, '-a', input_filename]
    if pbs:
        cmd.extend(['--sshloginfile $PBS_NODEFILE', '--workdir $PWD'])
    if isinstance(command, list):
        command = ' '.join(command)
    cmd.extend([command, '>', output_filename])
    return ' '.join(map(str, cmd))
