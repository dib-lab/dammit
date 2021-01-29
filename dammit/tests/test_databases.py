# Copyright (C) 2015-2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the BSD license.  See the LICENSE file for details.

import os
from os import path

from .utils import run
import pytest

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


@pytest.mark.skip
class TestDatabases():
    '''Tests for the dammit databases subcommand and the
    databases module. Assumes that DAMMIT_DB_DIR has been exported
    and a full install has already been performed.
    '''

    def test_default_database_dir_noenv(self):
        saved_env = os.environ.get('DAMMIT_DB_DIR', None)
        try:
            del os.environ['DAMMIT_DB_DIR']
        except KeyError:
            pass

        result = databases.default_database_dir(self.logger)
        expected = path.expanduser('~/.dammit/databases')
        assert result == expected

        if saved_env is not None:
            os.environ['DAMMIT_DB_DIR'] = saved_env

    def test_default_database_dir_env(self):
        saved_env = os.environ.get('DAMMIT_DB_DIR', None)
        try:
            del os.environ['DAMMIT_DB_DIR']
        except KeyError:
            pass

        expected = path.expanduser('~/.dammit/databases')
        os.environ['DAMMIT_DB_DIR'] = expected

        result = databases.default_database_dir(self.logger)

        assert result == expected

        del os.environ['DAMMIT_DB_DIR']
        if saved_env is not None:
            os.environ['DAMMIT_DB_DIR'] = saved_env

    @pytest.mark.requires_databases
    def test_dammit_databases_check(self):
        '''Test the database check subcommand.
        '''

        status, out, err = run(['databases'])
        assert 'All database tasks up-to-date.' in out

    def test_dammit_databases_check_fail(self, tmpdir):
        '''Test that the database check fails properly.
        '''
            
        args = ['databases', '--database-dir', str(tmpdir)]
        status, out, err = run(args, fail_ok=True)
        assert 'Must install databases' in out
        assert status == 2

    @pytest.mark.huge
    def test_dammit_database_install(self, tmpdir):
        '''Run a database installation (very long).
        '''
        args = ['databases', '--install', '--database-dir', str(tmpdir)]
        status, out, err = run(args)
