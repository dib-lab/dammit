# Copyright (C) 2015-2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the BSD license.  See the LICENSE file for details.

import os
import pytest

from ope.io import gff3

from .utils import run


def compare_gff(fn_a, fn_b):
    df_a = gff3.GFF3Parser(fn_a).read().sort_values(['seqid', 'start', 'end', 'ID', 'Target'])
    df_a.reset_index(inplace=True, drop=True)
    df_b = gff3.GFF3Parser(fn_b).read().sort_values(['seqid', 'start', 'end', 'ID', 'Target'])
    df_b.reset_index(inplace=True, drop=True)
    
    print('First DF:', df_a, '\n', '=' * 40)
    print('Second DF:', df_b, '\n', '=' * 40)
    return df_a.equals(df_b)


class TestDammitAnnotate:
    '''Integration: dammit run'''

    @pytest.mark.long
    @pytest.mark.requires_databases
    @pytest.mark.parametrize('n_threads', (1,4))
    def test_annotate_default(self, tmpdir, datadir, n_threads):
        '''--n-threads [N] annotate [INPUT.fa]
        '''
        with tmpdir.as_cwd():
            transcripts = datadir('pom.20.fa')
            exp_gff3 = datadir('pom.20.dammit.gff3')
            exp_fasta = datadir('pom.20.dammit.fasta')

            args = ['run', '--n-threads', str(n_threads), 'annotate', transcripts]
            status, out, err = run(*args)

            outdir = 'pom.20.dammit'
            gff3_fn = os.path.join(outdir, 'pom.20.dammit.gff3')
            fasta_fn = os.path.join(outdir, 'pom.20.dammit.fasta')

            print(os.listdir(outdir))
            print(gff3_fn, fasta_fn)
            assert compare_gff(gff3_fn, exp_gff3)
            assert open(fasta_fn).read() == open(exp_fasta).read()

    @pytest.mark.long
    @pytest.mark.requires_databases
    def test_annotate_evalue(self, tmpdir, datadir):
        '''annotate --global-evalue [EVALUE] [INPUT.fa]
        '''

        with tmpdir.as_cwd():
            transcripts = datadir('pom.20.fa')
            exp_gff3 = datadir('pom.20.dammit.evalue10.gff3')
            exp_fasta = datadir('pom.20.dammit.evalue10.fasta')

            args = ['run', 'annotate', transcripts, '--global-evalue', '10.0']
            status, out, err = run(*args)

            outdir = 'pom.20.dammit'
            gff3_fn = os.path.join(outdir, 'pom.20.dammit.gff3')
            fasta_fn = os.path.join(outdir, 'pom.20.dammit.fasta')

            assert compare_gff(gff3_fn, exp_gff3)
            assert open(fasta_fn).read() == open(exp_fasta).read()
  
    @pytest.mark.parametrize('n_threads', (1,4))
    def test_annotate_user_database(self, tmpdir, datadir, n_threads):
        '''--pipeline quick annotate --user-database [PEP.fa] [INPUT.fa]
        '''

        with tmpdir.as_cwd():
            transcripts = datadir('pom.20.fa')
            pep = datadir('pep.fa')
            exp_gff3 = datadir('pom.20.udb.dammit.gff3')
            exp_fasta = datadir('pom.20.udb.dammit.fasta')

            args = ['run', '--n-threads', str(n_threads), '--pipeline', 'quick', 'annotate',
                    transcripts, '--user-database', pep]
            status, out, err = run(*args)

            outdir = 'pom.20.dammit'
            gff3_fn = os.path.join(outdir, 'pom.20.dammit.gff3')
            fasta_fn = os.path.join(outdir, 'pom.20.dammit.fasta')

            assert status == 0
            assert compare_gff(gff3_fn, exp_gff3)
            assert open(fasta_fn).read() == open(exp_fasta).read()


#    def test_annotate_multiple_user_databases(self, tmpdir, datadir):
#        '''--pipeline quick annotate --user-database [PEP1.fa] --user-database [PEP2.fa] [INPUT.fa]
#        '''
#
#        with tmpdir.as_cwd():
#            transcripts = datadir('pom.single.fa')
#            pep = datadir('pep.fa')
#            pep2 = datadir('odb_subset.fa')
#            exp_gff3 = datadir('pom.single.fa.dammit.gff3.udb')
#            exp_fasta = datadir('pom.single.fa.dammit.fasta.udb')
#
#            args = ['annotate', '--quick',
#                    transcripts, '--user-databases', pep, pep2,
#                    '--verbosity', '2']
#            status, out, err = run(args)
#
#            outdir = '{0}.dammit'.format(transcripts)
#
#            assert status == 0

    def test_annotate_basename(self, tmpdir, datadir):
        '''--pipeline quick annotate --base-name [NAME] [INPUT.fa]
        '''

        with tmpdir.as_cwd():
            transcripts = datadir('pom.20.fa')

            args = ['run', '--pipeline', 'quick', 'annotate',
                    transcripts, '--base-name', 'Test']
            status, out, err = run(*args)
            assert status == 0

            fn = os.path.join('pom.20.dammit', 'pom.20.fasta')
            assert os.path.isfile(fn)

            contents = open(fn).read()
            assert 'Test_0' in contents


def test_annotate_dbdir_fail(tmpdir, datadir):
    '''Test annotation with a faulty database directory.
    '''

    with tmpdir.as_cwd():
        transcripts = datadir('pom.20.fa')

        args = ['run', '--database-dir', '.', 'annotate', transcripts]
        status, out, err = run(*args, fail_ok=True)
        print(status, out, err)
        assert 'you probably need to install the dammit databases' in err
        assert status == 1


def test_annotate_dbdir(tmpdir, datadir):
    '''Test that --database-dir works.
    '''

    with tmpdir.as_cwd():
        transcripts = datadir('pom.20.fa')

        database_dir = os.environ['DAMMIT_DB_DIR']
        args = ['run', '--database-dir', database_dir, '--pipeline', 'quick', 'annotate',  transcripts]
        status, out, err = run(*args)
        assert status == 0