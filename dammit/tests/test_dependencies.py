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


class TestDependencies(TestCase):

    @classmethod
    def setup_class(cls):
        pass
        #cls.PATH_BACKUP = os.environ['PATH']
        #print('Setup. PATH backup:', cls.PATH_BACKUP)

    @classmethod
    def teardown_class(cls):
        os.environ['PATH'] = PATH_BACKUP

    @staticmethod
    def add_execs_to_path(tempdir):
        os.environ['PATH'] = ''
        for exe in execs:
            fn = os.path.join(tempdir, exe)
            open(fn, 'a').close()
            os.chmod(fn, stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
        os.environ['PATH'] += os.pathsep + tempdir
        print('Added execs:', os.environ['PATH'])
        print('dir:', os.listdir(tempdir))

    def test_check_default_deps_nodeps(self):
        '''Test run_checks() when the dependencies aren't on the PATH.
        '''

        os.environ['PATH'] = ''
        handler = dependencies.DependencyHandler()
        for name, stat, msg in handler.run_checks():
            self.assertIn(name, names)
            self.assertFalse(stat)

    def test_check_system_path_alldeps(self):
        '''Test run_checks() when the dependencies are on the PATH.
        '''
        handler = dependencies.DependencyHandler()
        with TemporaryDirectory() as tempdir:
            TestDependencies.add_execs_to_path(tempdir)

            for name, stat, msg in handler.run_checks():
                self.assertIn(name, names)
                if name != 'LAST':
                    self.assertTrue(stat, msg=name + ' ' + msg)
    
    def test_handle_nodeps(self):
        os.environ['PATH'] = ''
        handler = dependencies.DependencyHandler()
        missing = handler.handle()
        for name in missing:
            self.assertIn(name, names)
        self.assertEquals(set(missing) - set(names), set())

    def test_do_check_alldeps(self):
        handler = dependencies.DependencyHandler()
        with TemporaryDirectory() as tempdir:
            TestDependencies.add_execs_to_path(tempdir)
            missing = handler.handle()
            self.assertEquals(missing, ['LAST'])

class TestDammitDependencies(TestCase):

    def setUp(self):
        from itertools import count
        gff.ID_GEN = count()


    def test_dammit_dependencies(self):
        '''Test the dependencies subcommand.
        '''

        status, out, err = run(['dependencies'])
        self.assertEquals(status, 0)


