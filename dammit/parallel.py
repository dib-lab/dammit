#!/usr/bin/env python
import os
from .utils import which, doit_task


def parallel_fasta(input_filename, output_filename, command, n_jobs, pbs=False):
    '''Given an input FASTA source, target, shell command, and number of jobs,
    construct a gnu-parallel command to act on the sequences.

    Args:
        input_filename (str): The source FASTA.
        output_filename (str): The target.
        command (list): The shell command (in subprocess format).
        n_jobs (int): Number of cores or nodes to split to.
        pbs (bool): If True, add the appropriate flags for running on a
        multinode system. Note that this means the user needs to export
        $PBS_NODEFILE before running dammit.
    Returns:
        str: The constructed shell command.

    '''

    exc = which('parallel')
    cmd = ['cat', input_filename, '|', exc, '--progress', '--pipe', '-L', 2, '-N', 400,
           '--gnu', '-j', n_jobs, '-a', input_filename]
    if pbs:
        cmd.extend(['--sshloginfile $PBS_NODEFILE', '--workdir $PWD'])
    if isinstance(command, list):
        command = ' '.join(command)
    cmd.extend([command, '>', output_filename])
    return ' '.join(map(str, cmd))
