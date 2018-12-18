# Copyright (C) 2015-2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the BSD license.  See the LICENSE file for details.

import os
import stat
import pandas as pd
import pytest
import traceback
import sys

from dammit.app import DammitApp
from dammit.fileio import gff3

from utils import datadir, runscript

PATH_BACKUP = os.environ['PATH']

def run(args, **kwargs):
    return runscript('dammit', args, **kwargs)


def compare_gff(fn_a, fn_b):
    df_a = gff3.GFF3Parser(fn_a).read().sort_values('ID')
    df_a.reset_index(inplace=True, drop=True)
    df_b = gff3.GFF3Parser(fn_b).read().sort_values('ID')
    df_b.reset_index(inplace=True, drop=True)
    
    print('First DF:', df_a, '\n', '=' * 40)
    print('Second DF:', df_b, '\n', '=' * 40)
    return df_a.equals(df_b)


class TestDammitAnnotate:

    def setup_method(self):
        '''This was fun to diagnose! Because the acceptance tests are actually
        being executed in the current namespace with eval, this global generator
        was retaining its state between tests. Oops!
        '''

        gff3.next_ID = gff3.id_gen_wrapper()
        self.maxDiff = None

    @pytest.mark.long
    @pytest.mark.requires_databases
    def test_annotate_basic(self, tmpdir, datadir):
        '''Run a basic annotation and verify the results.
        '''

        with tmpdir.as_cwd():
            transcripts = datadir('pom.single.fa')
            exp_gff3 = datadir('pom.single.fa.dammit.gff3')
            exp_fasta = datadir('pom.single.fa.dammit.fasta')

            args = ['annotate', transcripts]
            status, out, err = run(args)

            outdir = '{0}.dammit'.format(transcripts)
            gff3_fn = os.path.join(outdir, 'pom.single.fa.dammit.gff3')
            fasta_fn = os.path.join(outdir, 'pom.single.fa.dammit.fasta')

            print(os.listdir(outdir))
            print(gff3_fn, fasta_fn)
            assert compare_gff(gff3_fn, exp_gff3)
            assert open(fasta_fn).read() == open(exp_fasta).read()

    @pytest.mark.long
    @pytest.mark.requires_databases
    def test_annotate_full(self, tmpdir, datadir):
        '''Run a full annotation and verify the results.
        '''

        with tmpdir.as_cwd():
            transcripts = datadir('pom.single.fa')
            exp_gff3 = datadir('pom.single.fa.dammit.gff3.full')
            exp_fasta = datadir('pom.single.fa.dammit.fasta.full')

            args = ['annotate', transcripts, '--full']
            status, out, err = run(args)

            outdir = '{0}.dammit'.format(transcripts)
            gff3_fn = os.path.join(outdir, 'pom.single.fa.dammit.gff3')
            fasta_fn = os.path.join(outdir, 'pom.single.fa.dammit.fasta')

            assert compare_gff(gff3_fn, exp_gff3)
            assert open(fasta_fn).read() == open(exp_fasta).read()

    @pytest.mark.long
    @pytest.mark.requires_databases
    def test_annotate_threaded(self, tmpdir, datadir):
        '''Test the --n_threads argument.
        '''

        with tmpdir.as_cwd():
            transcripts = datadir('pom.single.fa')
            args = ['annotate', transcripts, '--n_threads', '2']
            status, out, err = run(args)

    @pytest.mark.long
    @pytest.mark.requires_databases
    def test_annotate_evalue(self, tmpdir, datadir):
        '''Test the --evalue argument and verify the results.
        '''

        with tmpdir.as_cwd():
            transcripts = datadir('pom.single.fa')
            exp_gff3 = datadir('pom.single.fa.dammit.gff3.evalue10')
            exp_fasta = datadir('pom.single.fa.dammit.fasta.evalue10')

            args = ['annotate', transcripts, '--evalue', '10.0']
            status, out, err = run(args)

            outdir = '{0}.dammit'.format(transcripts)
            gff3_fn = os.path.join(outdir, 'pom.single.fa.dammit.gff3')
            fasta_fn = os.path.join(outdir, 'pom.single.fa.dammit.fasta')

            assert compare_gff(gff3_fn, exp_gff3)
            assert open(fasta_fn).read() == open(exp_fasta).read()
  

    def test_annotate_outdir(self, tmpdir, datadir):
        '''Test that the --output-dir argument works.
        '''

        with tmpdir.as_cwd():
            transcripts = datadir('pom.single.fa')
            outdir = 'test_out'
            args = ['annotate', '--quick', transcripts, '-o', outdir]
            status, out, err = run(args)
            fn = os.path.join(outdir, os.path.basename(transcripts))
            assert os.path.isfile(fn)

    def test_annotate_dbdir_fail(self, tmpdir, datadir):
        '''Test annotation with a faulty database directory.
        '''

        with tmpdir.as_cwd():
            transcripts = datadir('pom.single.fa')

            args = ['annotate', transcripts, '--database-dir', '.']
            status, out, err = run(args, fail_ok=True)
            assert 'install databases to continue' in out
            assert status == 2

    def test_annotate_dbdir(self, tmpdir, datadir):
        '''Test that --database-dir works.
        '''

        with tmpdir.as_cwd():
            transcripts = datadir('pom.single.fa')

            db_dir = os.environ['DAMMIT_DB_DIR']
            args = ['annotate', '--quick', transcripts, '--database-dir', db_dir]
            status, out, err = run(args)

    def test_annotate_user_databases(self, tmpdir, datadir):
        '''Test that a user database works.
        '''

        with tmpdir.as_cwd():
            transcripts = datadir('pom.single.fa')
            pep = datadir('pep.fa')
            exp_gff3 = datadir('pom.single.fa.dammit.gff3.udb')
            exp_fasta = datadir('pom.single.fa.dammit.fasta.udb')

            args = ['annotate', '--quick',
                    transcripts, '--user-databases', pep,
                    '--verbosity', '2']
            status, out, err = run(args)

            outdir = '{0}.dammit'.format(transcripts)
            gff3_fn = os.path.join(outdir, 'pom.single.fa.dammit.gff3')
            fasta_fn = os.path.join(outdir, 'pom.single.fa.dammit.fasta')

            assert status == 0
            assert compare_gff(gff3_fn, exp_gff3)
            assert open(fasta_fn).read() == open(exp_fasta).read()

    def test_annotate_multiple_user_databases(self, tmpdir, datadir):
        '''Test that multiple user databases work.
        '''

        with tmpdir.as_cwd():
            transcripts = datadir('pom.single.fa')
            pep = datadir('pep.fa')
            pep2 = datadir('odb_subset.fa')
            exp_gff3 = datadir('pom.single.fa.dammit.gff3.udb')
            exp_fasta = datadir('pom.single.fa.dammit.fasta.udb')

            args = ['annotate', '--quick',
                    transcripts, '--user-databases', pep, pep2,
                    '--verbosity', '2']
            status, out, err = run(args)

            outdir = '{0}.dammit'.format(transcripts)

            assert status == 0

    def test_annotate_name(self, tmpdir, datadir):
        '''Test the --name argument.
        '''

        with tmpdir.as_cwd():
            transcripts = datadir('pom.single.fa')

            args = ['annotate', '--quick',
                    transcripts, '--name', 'Test']
            status, out, err = run(args)

            outdir = '{0}.dammit'.format(transcripts)
            fn = os.path.join(outdir, os.path.basename(transcripts))
            assert os.path.isfile(fn)

            contents = open(fn).read()
            assert 'Test_0' in contents

            assert status == 0
    
    
    def test_annotate_no_rename(self, tmpdir, datadir):
        '''Test the --no-rename argument.
        '''

        with tmpdir.as_cwd():
            transcripts = datadir('pom.single.fa')

            args = ['annotate', transcripts, '--no-rename']
            status, out, err = run(args)

            outdir = '{0}.dammit'.format(transcripts)
            fn = os.path.join(outdir, os.path.basename(transcripts))
            assert os.path.isfile(fn)

            contents = open(fn).read()
            assert 'SPAC212' in contents

            assert status == 0

