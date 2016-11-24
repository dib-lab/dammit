#!/usr/bin/env python
from __future__ import print_function

import subprocess

from .utils import DependentTask, InstallationError
from ..utils import which, doit_task
from ..profile import profile_task

from shmlast.last import lastdb_task
from shmlast.last import lastal_task


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
        return lastdb_task(*args, **kwargs)


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
        return lastal_task(*args, **kwargs)
