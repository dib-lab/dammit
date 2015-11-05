from unittest import TestCase

import os
import stat

from dammit.app import DammitApp
from dammit import dependencies

from utils import TemporaryDirectory

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

class TestDammitApp(TestCase):

    def test_main_module(self):
        pass
        

class TestDependencies(TestCase):

    @classmethod
    def setup_class(cls):
        cls.PATH_BACKUP = os.environ['PATH']
        print('Setup. PATH backup:', cls.PATH_BACKUP)

    def teardown(self):
        '''Most of the functions need to modify PATH to test properly.
        Make sure to fix it when we're done.
        '''

        os.environ['PATH'] = self.PATH_BACKUP

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

    def test_check_system_path_nodeps(self):
        '''Test check_system_path when the dependencies aren't on the PATH.
        '''

        os.environ['PATH'] = ''
        handler = dependencies.DependencyHandler()
        deps = handler.check_system_path()
        for name, stat in deps.iteritems():
            self.assertIn(name, names)
            self.assertFalse(stat)

    def test_check_system_path_alldeps(self):
        '''Test check_system_path when the dependencies are on the PATH.
        '''
        handler = dependencies.DependencyHandler()
        with TemporaryDirectory() as tempdir:
            TestDependencies.add_execs_to_path(tempdir)

            deps = handler.check_system_path()
            for name, stat in deps.iteritems():
                self.assertIn(name, names)
                self.assertTrue(stat, msg=name)
    
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
            self.assertEquals(missing, [])

class TestDatabases(TestCase):

    @classmethod
    def setup_class(cls):
        cls.db_dir = 'test_db_dir'


            
