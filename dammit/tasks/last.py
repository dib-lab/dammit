# Copyright (C) 2015-2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the BSD license.  See the LICENSE file for details.

import subprocess

from dammit.tasks.utils import DependentTask, InstallationError
from dammit.utils import which, doit_task
from dammit.profile import profile_task

from shmlast import last


def check_version(exc):
    try:
        version_string = subprocess.check_output([exc, '--version'])
    except subprocess.CalledProcessError as e:
        raise InstallationError('Error checking {0} version: [{1}] {2}'.format(exc,
                                                                     e.returncode, 
                                                                     e.output))
    except OSError as e:
        raise InstallationError('Error checking {0} version: [{1}] {2}'.format(exc,
                                                                     e.errno, 
                                                                     str(e)))
    else:
        _, version = version_string.strip().split()
        return int(version)


class LastDBTask(DependentTask):

    def deps(self):
        lastdb = which('lastdb')
        if lastdb is None:
            raise InstallationError('lastdb not found.')
        version = check_version(lastdb)
        if version < 600:
            raise InstallationError('lastdb version {0} < 600, please'\
                                    ' update'.format(version))
        if self.logger:
            self.logger.debug('lastdb:' + lastdb)

        return lastdb

    def task(self, *args, **kwargs):
        exc = self.deps()
        return last.lastdb_task(*args, **kwargs)


class LastalTask(DependentTask):

    def deps(self):
        lastal = which('lastdb')
        if lastal is None:
            raise InstallationError('lastal not found.')
        version = check_version(lastal)
        if version < 600:
            raise InstallationError('lastal version {0} < 600, please'\
                                    ' update'.format(version))
        if self.logger:
            self.logger.debug('lastal:' + lastal)

        return lastal
    
    def task(self, *args, **kwargs):
        exc = self.deps()
        return last.lastal_task(*args, **kwargs)
