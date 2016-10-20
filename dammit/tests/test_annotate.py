from unittest import TestCase

from nose.plugins.attrib import attr
import os
import stat
import pandas as pd

from dammit.app import DammitApp
from dammit import dependencies
from dammit.fileio import gff3

from utils import TemporaryDirectory, TestData, runscript

PATH_BACKUP = os.environ['PATH']

def run(args, **kwargs):
    return runscript('dammit', args, **kwargs)


def compare_gff(fn_a, fn_b):
    df_a = gff3.GFF3Parser(fn_a).read().sort_values('ID')
    df_a.reset_index(inplace=True, drop=True)
    df_b = gff3.GFF3Parser(fn_b).read().sort_values('ID')
    df_b.reset_index(inplace=True, drop=True)
    
    print(df_a)
    print(df_b)
    return df_a.equals(df_b)


@attr('long')
class TestDammitAnnotate(TestCase):

    def setUp(self):
        '''This was fun to diagnose! Because the acceptance tests are actually
        being executed in the current namespace with eval, this global generator
        was retaining its state between tests. Oops!
        '''

        gff3.next_ID = gff3.id_gen_wrapper()
        self.maxDiff = None

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

            print(os.listdir(td))
            print(os.listdir(outdir))
            print(gff3_fn, fasta_fn)
            assert compare_gff(gff3_fn, exp_gff3)
            self.assertEquals(open(fasta_fn).read(), open(exp_fasta).read())

    def test_annotate_full(self):
        '''Run a full annotation and verify the results.
        '''

        with TemporaryDirectory() as td,\
             TestData('pom.single.fa', td) as transcripts,\
             TestData('pom.single.fa.dammit.gff3.full', td) as exp_gff3,\
             TestData('pom.single.fa.dammit.fasta.full', td) as exp_fasta:

            args = ['annotate', transcripts, '--full']
            status, out, err = run(args, in_directory=td)

            outdir = '{0}.dammit'.format(transcripts)
            gff3_fn = os.path.join(outdir, 'pom.single.fa.dammit.gff3')
            fasta_fn = os.path.join(outdir, 'pom.single.fa.dammit.fasta')

            assert compare_gff(gff3_fn, exp_gff3)
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

            assert compare_gff(gff3_fn, exp_gff3)
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
            self.assertIn('install databases to continue', out)
            self.assertEquals(status, 2)

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
             TestData('pep.fa', td) as pep, \
             TestData('pom.single.fa.dammit.gff3.udb', td) as exp_gff3,\
             TestData('pom.single.fa.dammit.fasta.udb', td) as exp_fasta:

            args = ['annotate', transcripts, '--user-databases', pep,
                    '--verbosity', '2']
            status, out, err = run(args, in_directory=td)

            outdir = '{0}.dammit'.format(transcripts)
            gff3_fn = os.path.join(outdir, 'pom.single.fa.dammit.gff3')
            fasta_fn = os.path.join(outdir, 'pom.single.fa.dammit.fasta')

            assert compare_gff(gff3_fn, exp_gff3)
            self.assertEquals(open(fasta_fn).read(), open(exp_fasta).read())

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


