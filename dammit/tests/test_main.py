from unittest import TestCase

from nose.plugins.attrib import attr
import os
import stat

from dammit.app import DammitApp
from dammit import dependencies


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


@attr('acceptance')
class TestDammit(TestCase):

    def test_dammit_version(self):
        '''Test the dammit --version command.
        '''

        from dammit import __version__
        status, out, err = run(['--version'])
        self.assertEquals(status, 0)
        self.assertEquals(err.strip(), 'dammit {0}'.format(__version__))

    def test_dammit_dependencies(self):
        '''Test the dammit dependencies subcommand.
        '''

        status, out, err = run(['dependencies'])
        self.assertEquals(status, 0)

    def test_dammit_databases_check(self):
        '''Test dammit database check subcommand.
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

    def test_annotate_basic(self):
        with TemporaryDirectory() as td,\
             TestData('pom.single.fa', td) as transcripts:

            args = ['annotate', transcripts]
            status, out, err = run(args, in_directory=td)

    def test_annotate_threaded(self):
        with TemporaryDirectory() as td,\
             TestData('pom.single.fa', td) as transcripts:

            args = ['annotate', transcripts, '--n_threads', '2']
            status, out, err = run(args, in_directory=td)


    def test_annotate_evalue(self):
        with TemporaryDirectory() as td,\
             TestData('pom.single.fa', td) as transcripts:

            args = ['annotate', transcripts, '--evalue', '.01']
            status, out, err = run(args, in_directory=td)
  

    def test_annotate_outdir(self):
        with TemporaryDirectory() as td,\
             TestData('pom.single.fa', td) as transcripts:

            outdir = os.path.join(td, 'test_out')
            args = ['annotate', transcripts, '-o', outdir]
            status, out, err = run(args)


class TestDammitApp(TestCase):

    def test_main_module(self):
        pass
        

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


            
