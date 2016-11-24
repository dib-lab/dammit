#!/usr/bin/env python
from __future__ import print_function

import logging
import os
import subprocess
import sys

from doit.dependency import Dependency, SqliteDB

from . import ui
from .utils import which
from .tasks.shell import get_download_and_untar_task


class DependencyHandler(object):

    def __init__(self):
        self.logger = logging.getLogger(self.__class__.__name__)
        self.checks = {}

    def register_dependency_check(self, name, function):
        '''Add a new dependency to the handler.

        Args:
            name (str): The name of the dependency.
            function (callable): A function which returns a tuple
                of (boolean, string) for the status and a message, and
                takes a Logger object.
        '''
        self.logger.debug('Add check function for {0}'.format(name))
        self.checks[name] = function

    def clear_dependency_checks(self):
        self.logger.debug('Clearing {0} dependency checks from'\
                          ' handler'.format(len(self.checks)))
        self.checks = {}

    def get_status(self, name):
        self.logger.debug('Run dependency function on {0}'.format(name))
        try:
            check_function = self.checks[name]
        except KeyError:
            self.logger.error('Check not found: {0}'.format(name))
        status, msg = check_function(self.logger)
        self.logger.debug('{0} has status {1}, msg: {2}'.format(name, status, msg))

        return status, msg

    def get_all_statuses(self):
        fulfilled = True
        unfulfilled = {}
        for name, _ in self.checks.items():
            status, msg = self.get_status(name)
            if status is False:
                fulfilled = False
                unfulfilled[name] = msg
        return fulfilled, unfulfilled

    def print_all_statuses(self, out=sys.stdout):
        is_fulfilled, unfulfilled = self.get_all_statuses()
        if is_fulfilled:
            print(ui.paragraph('*All dependencies fulfilled.*'), file=out)
        else:
            print('\nSome dependencies unfulfilled:', file=out)
            print(ui.listing(unfulfilled), file=out)
        return is_fulfilled, unfulfilled
            
    def check_or_fail(self, out=sys.stdout):
        print(ui.header('Dependency Check', level=3), file=out)
        is_fulfilled, unfulfilled = self.print_all_statuses(out=out)
        if not is_fulfilled:
            print(ui.paragraph('Must install dependencies to continue.'\
                              ' to do so, follow the directions in the'\
                              ' documentation at '\
                              'http://www.camillescott.org/dammit/installing.html'),
                 file=out)
            sys.exit(1)


def get_handler():
    return register_builtin_checks(DependencyHandler())


def register_builtin_checks(handler):
    checks = {'TransDecoder': check_transdecoder,
               'LAST': check_last}
    for name, func in checks.items():
        handler.register_dependency_check(name, func)
    return handler


def check_transdecoder(logger):
    longorfs = which('TransDecoder.LongOrfs')
    predict = which('TransDecoder.Predict')
    if longorfs is None or predict is None:
        return False, 'Not found on $PATH'
    else:
        logger.debug('TransDecoder.LongOrfs:' + longorfs)
        logger.debug('TransDecoder.Predict:' + predict)
        return True, os.path.dirname(longorfs)


def check_last(logger):
    lastdb = which('lastdb')
    lastal = which('lastal')
    if lastdb is None or lastal is None:
        return False, 'Not found on $PATH'
    else:
        try:
            version_string = subprocess.check_output(['lastal', '--version'])
        except subprocess.CalledProcessError as e:
            return False, 'Error checking lastal version: [{0}] {1}'.format(e.returncode, e.output)
        except OSError as e:
            return False, 'Error checking lastal version: [{0}] {1}'.format(e.errno, str(e))
        else:
            _, version = version_string.strip().split()
            if int(version) < 600:
                return False, 'lastal version {0} < 600, please update'.format(version)
            logger.debug('lastal:' + lastal)
            logger.debug('latsdb:' + lastdb)
            return True, os.path.dirname(lastal)


