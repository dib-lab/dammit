import os
from os import path
import logging
import stat
import pandas as pd

from dammit import databases
from dammit.meta import get_config

from utils import datadir, runscript, logger
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

PATH_BACKUP = os.environ['PATH']

def run(args, **kwargs):
    return runscript('dammit', args, **kwargs)


class TestDatabases():
    '''Tests for the dammit databases subcommand and the
    databases module. Assumes that DAMMIT_DB_DIR has been exported
    and a full install has already been performed.
    '''

    def setup_method(self):
        self.logger = logging.getLogger('tests.test_databases.TestDatabases')
        class Args(object):
            pass
        self.args = Args()
        self.args.database_dir = databases.default_database_dir(self.logger)
        self.args.verbosity = 2
        self.args.full = False
        self.args.n_threads = 1
        self.args.busco_group = 'metazoa'

        self.config, self.databases = get_config()
        self.config.update(vars(self.args))
        self.handler = databases.get_handler(self.config)
        databases.build_default_pipeline(self.handler,
                                         self.config,
                                         self.databases,
                                         with_uniref=False)


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

    @pytest.mark.requires_databases
    def test_check_or_fail_succeed(self):
        '''Check that check_or_fail succeeds properly.
        '''
        try:
            databases.check_or_fail(self.handler)
        except SystemExit:
            assert False, 'Should not have exited'

    def test_check_or_fail_fail(self, tmpdir):
        config = self.config.copy()
        config['database_dir'] = str(tmpdir)
        
        handler = databases.get_handler(config)
        databases.build_default_pipeline(handler,
                                         config,
                                         self.databases,
                                         with_uniref=False)
        with pytest.raises(SystemExit):
            databases.check_or_fail(handler)
