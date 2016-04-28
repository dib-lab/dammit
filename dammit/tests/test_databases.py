from unittest import TestCase

from nose.plugins.attrib import attr
import os
import stat
import pandas as pd

from dammit import databases
from dammit.meta import get_config

from utils import TemporaryDirectory, TestData, runscript

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
        class Args(object):
            pass
        self.args = Args()
        self.args.database_dir = None
        self.args.verbosity = 2
        self.args.full = False
        self.config, self.databases = get_config()
        self.handler = databases.get_handler(self.args, self.config,
                                             self.databases)

    def test_dammit_databases_check(self):
        '''Test the database check subcommand.
        '''

        status, out, err = run(['databases'])
        self.assertIn('All databases prepared!', err)

    def test_dammit_databases_check_fail(self):
        '''Test that the database check fails properly.
        '''
        with TemporaryDirectory() as td:
            
            args = ['databases', '--database-dir', td]
            status, out, err = run(args, fail_ok=True)
            self.assertIn('prep incomplete', err)
            self.assertEquals(status, 1)

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
