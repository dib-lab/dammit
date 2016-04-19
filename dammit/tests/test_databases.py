from unittest import TestCase

from nose.plugins.attrib import attr
import os
import stat
import pandas as pd

from dammit.app import DammitApp
from dammit import dependencies
from dammit import gff

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


class TestDammitDatabases(TestCase):

    def setUp(self):
        from itertools import count
        gff.ID_GEN = count()

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
    def test_dammit_database_install_full(self):
        '''Run a full database installation (very long).
        '''

        with TemporaryDirectory() as td:
            args = ['databases', '--install', '--database-dir', td]
            status, out, err = run(args)


