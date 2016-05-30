from unittest import TestCase

from nose.plugins.attrib import attr
import os
from os import path
import logging
import stat
import pandas as pd

from dammit import databases
from dammit.meta import get_config

from utils import TemporaryDirectory, TestData, runscript, logger

names = ['TransDecoder',
         'crb-blast',
         'LAST',
         'BUSCO',
         'BLAST+',
         'Infernal',
         'HMMER']

execs = ['hmmscan',
         'hmmpress',
         'cmscan',
         'cmpress',
         'blastp',
         'blastx',
         'tblastn',
         'makeblastdb',
         'BUSCO_v1.1b1.py',
         'TransDecoder.LongOrfs',
         'TransDecoder.Predict',
         'lastal',
         'lastdb',
         'crb-blast']

PATH_BACKUP = os.environ['PATH']

def run(args, **kwargs):
    return runscript('dammit', args, **kwargs)


class TestDatabases(TestCase):
    '''Tests for the dammit databases subcommand and the
    databases module. Assumes that DAMMIT_DB_DIR has been exported
    and a full install has already been performed.
    '''

    def setUp(self):
        self.logger = logging.getLogger('TestCase.TestDatabases')
        class Args(object):
            pass
        self.args = Args()
        self.args.database_dir = databases.default_database_dir(self.logger)
        self.args.verbosity = 2
        self.args.full = False
        self.config, self.databases = get_config()
        self.handler = databases.get_handler(self.args, self.config,
                                             self.databases)


    def test_default_database_dir_noenv(self):
        saved_env = os.environ.get('DAMMIT_DB_DIR', None)
        try:
            del os.environ['DAMMIT_DB_DIR']
        except KeyError:
            pass

        result = databases.default_database_dir(self.logger)
        expected = path.expanduser('~/.dammit/databases')
        self.assertTrue(result == expected)

        if saved_env is not None:
            os.environ['DAMMIT_DB_DIR'] = saved_env

    def test_default_database_dir_env(self):
        expected = path.expanduser('~/.dammit/databases')
        os.environ['DAMMIT_DB_DIR'] = expected

        result = databases.default_database_dir(self.logger)

        self.assertTrue(result == expected)

        del os.environ['DAMMIT_DB_DIR']

    def test_dammit_databases_check(self):
        '''Test the database check subcommand.
        '''

        status, out, err = run(['databases'])
        self.assertIn('All tasks up-to-date!', out)

    def test_dammit_databases_check_fail(self):
        '''Test that the database check fails properly.
        '''
        with TemporaryDirectory() as td:
            
            args = ['databases', '--database-dir', td]
            status, out, err = run(args, fail_ok=True)
            self.assertIn('Out-of-date tasks', out)
            self.assertEquals(status, 2)

    @attr('huge')
    def test_dammit_database_install(self):
        '''Run a database installation (very long).
        '''
        with TemporaryDirectory() as td:
            args = ['databases', '--install', '--database-dir', td]
            status, out, err = run(args)

    def test_check_or_fail_succeed(self):
        '''Check that check_or_fail succeeds properly.
        '''
        try:
            databases.check_or_fail(self.handler)
        except SystemExit:
            assert False, 'Should not have exited'

    def test_check_or_fail_fail(self):
        with TemporaryDirectory() as td:
            self.args.database_dir = td
            handler = databases.get_handler(self.args, self.config,
                                                 self.databases)
            with self.assertRaises(SystemExit):
                databases.check_or_fail(handler)
