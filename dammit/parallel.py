#!/usr/bin/env python
import os
import subprocess

from .utils import which, doit_task
from .tasks.utils import InstallationError


def check_parallel(logger=None):
    parallel = which('parallel')
    if parallel is None:
        raise InstallationError('parallel not found.')
    else:
        try:
            version_string = subprocess.check_output(['parallel', '--version'])
        except subprocess.CalledProcessError as e:
            raise InstallationError('Error checking parallel '\
                                    'version: [{0}] {1}'.format(e.returncode, e.output))
        except OSError as e:
            raise InstallationError('Error checking parallel version: '\
                                    '[{0}] {1}'.format(e.errno, str(e)))
        else:
            version = version_string.strip().split()[2]
            if logger:
                logger.debug('parallel version:{0}'.format(version))
            if int(version) < 20150000:
                raise InstallationError('parallel version {0} < 20150000, '\
                                        'please update'.format(version))
            if logger:
                logger.debug('parallel:' + parallel)
            return parallel


def parallel_fasta(input_filename, output_filename, command, n_jobs,
                   pbs=False, check_dep=True, logger=None):
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
    
    exc = which('parallel') if not check_dep else check_parallel(logger=logger)
    cmd = ['cat', input_filename, '|', exc, '--progress', '--pipe', '-L', 2, '-N', 400,
           '--gnu', '-j', n_jobs, '-a', input_filename]
    if pbs:
        cmd.extend(['--sshloginfile $PBS_NODEFILE', '--workdir $PWD'])
    if isinstance(command, list):
        command = ' '.join(command)
    cmd.extend([command, '>', output_filename])
    return ' '.join(map(str, cmd))
