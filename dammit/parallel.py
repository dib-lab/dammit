# Copyright (C) 2015-2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the BSD license.  See the LICENSE file for details.

import os
import subprocess

from dammit.utils import which, doit_task
from dammit.tasks.utils import InstallationError


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
                   sshloginfile=None, check_dep=True, logger=None):
    '''Given an input FASTA source, target, shell command, and number of jobs,
    construct a gnu-parallel command to act on the sequences.

    Args:
        input_filename (str): The source FASTA.
        output_filename (str): The target.
        command (list): The shell command (in subprocess format).
        n_jobs (int): Number of cores or nodes to split to.
        sshloginfile (str): Path to file with node addresses.
        check_dep (bool): If True, check for the gnu-parallel executable.
        logger (logging.Logger): A logger to use.
    Returns:
        str: The constructed shell command.
    '''

    exc = which('parallel') if not check_dep else check_parallel(logger=logger)
    cmd = ['cat', input_filename, '|', exc, '--round-robin', '--pipe', '-L', 2,
           '-N', 10000, '--gnu']
    if sshloginfile is not None:
        cmd.extend(['--sshloginfile', sshloginfile, '--workdir $PWD'])
    else:
        cmd.extend(['-j', n_jobs])
    cmd.extend(['-a', input_filename])

    if isinstance(command, list):
        command = ' '.join(command)
    cmd.extend([command, '>', output_filename])
    return ' '.join(map(str, cmd))
