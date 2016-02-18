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


@attr('acceptance')
class TestDammit(TestCase):

    def setUp(self):
        '''This was fun to diagnose! Because the acceptance tests are actually
        being executed in the current namespace with eval, this global generator
        was retaining its state between tests. Oops!
        '''

        from itertools import count
        gff.ID_GEN = count()

    def test_dammit_version(self):
        '''Test the dammit --version command.
        '''

        from dammit import __version__
        status, out, err = run(['--version'])
        self.assertEquals(status, 0)
        self.assertEquals(err.strip(), 'dammit {0}'.format(__version__))

    def test_dammit_dependencies(self):
        '''Test the dependencies subcommand.
        '''

        status, out, err = run(['dependencies'])
        self.assertEquals(status, 0)

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

    @attr('long')
    def test_dammit_database_install_full(self):
        '''Run a full database installation (very long).
        '''

        with TemporaryDirectory() as td:
            args = ['databases', '--install', '--database-dir', td]
            status, out, err = run(args)

    def test_annotate_basic(self):
        '''Run a basic annotation and verify the results.
        '''

        with TemporaryDirectory() as td,\
             TestData('pom.single.fa', td) as transcripts,\
             TestData('pom.single.fa.dammit.gff3', td) as exp_gff3,\
             TestData('pom.single.fa.dammit.fasta', td) as exp_fasta:

            args = ['annotate', transcripts]
            status, out, err = run(args, in_directory=td)

            outdir = '{0}.dammit'.format(transcripts)
            gff3_fn = os.path.join(outdir, 'pom.single.fa.dammit.gff3')
            fasta_fn = os.path.join(outdir, 'pom.single.fa.dammit.fasta')

            self.assertEquals(open(gff3_fn).read(), open(exp_gff3).read())
            self.assertEquals(open(fasta_fn).read(), open(exp_fasta).read())

    def test_annotate_threaded(self):
        '''Test the --n_threads argument.
        '''

        with TemporaryDirectory() as td,\
             TestData('pom.single.fa', td) as transcripts:

            args = ['annotate', transcripts, '--n_threads', '2']
            status, out, err = run(args, in_directory=td)


    def test_annotate_evalue(self):
        '''Test the --evalue argument and verify the results.
        '''

        with TemporaryDirectory() as td,\
             TestData('pom.single.fa', td) as transcripts,\
             TestData('pom.single.fa.dammit.gff3.evalue10', td) as exp_gff3,\
             TestData('pom.single.fa.dammit.fasta.evalue10', td) as exp_fasta:

            args = ['annotate', transcripts, '--evalue', '10.0']
            status, out, err = run(args, in_directory=td)

            outdir = '{0}.dammit'.format(transcripts)
            gff3_fn = os.path.join(outdir, 'pom.single.fa.dammit.gff3')
            fasta_fn = os.path.join(outdir, 'pom.single.fa.dammit.fasta')

            self.assertEquals(open(gff3_fn).read(), open(exp_gff3).read())
            self.assertEquals(open(fasta_fn).read(), open(exp_fasta).read())
  

    def test_annotate_outdir(self):
        '''Test that the --output-dir argument works.
        '''

        with TemporaryDirectory() as td,\
             TestData('pom.single.fa', td) as transcripts:

            outdir = os.path.join(td, 'test_out')
            args = ['annotate', transcripts, '-o', outdir]
            status, out, err = run(args)
            fn = os.path.join(outdir, os.path.basename(transcripts))
            self.assertTrue(os.path.isfile(fn))

    def test_annotate_dbdir_fail(self):
        '''Test annotation with a faulty database directory.
        '''

        with TemporaryDirectory() as td, \
             TestData('pom.single.fa', td) as transcripts:

            args = ['annotate', transcripts, '--database-dir', td]
            status, out, err = run(args, in_directory=td, fail_ok=True)
            self.assertIn('Install databases to continue', err)
            self.assertEquals(status, 1)

    def test_annotate_dbdir(self):
        '''Test that --database-dir works.
        '''

        with TemporaryDirectory() as td,\
             TestData('pom.single.fa', td) as transcripts:

            db_dir = os.environ['DAMMIT_DB_DIR']
            args = ['annotate', transcripts, '--database-dir', db_dir]
            status, out, err = run(args, in_directory=td)

    def test_annotate_user_databases(self):
        '''Test that a user database works.
        '''

        with TemporaryDirectory() as td,\
             TestData('pom.single.fa', td) as transcripts,\
             TestData('pep.fa', td) as pep:

            args = ['annotate', transcripts, '--user-databases', pep]
            status, out, err = run(args, in_directory=td)

    def test_annotate_name(self):
        '''Test the --name argument.
        '''

        with TemporaryDirectory() as td,\
             TestData('pom.single.fa', td) as transcripts:

            args = ['annotate', transcripts, '--name', 'Test']
            status, out, err = run(args, in_directory=td)
            outdir = '{0}.dammit'.format(transcripts)
            fn = os.path.join(outdir, os.path.basename(transcripts))
            self.assertTrue(os.path.isfile(fn))
            contents = open(fn).read()
            self.assertIn('Test_0', contents)


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

class TestDatabases(TestCase):

    @classmethod
    def setup_class(cls):
        cls.db_dir = 'test_db_dir'


            
