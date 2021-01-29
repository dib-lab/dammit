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
    def test_default(self, tmpdir, datadir, n_threads):
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
    def test_evalue(self, tmpdir, datadir):
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
    def test_user_database(self, tmpdir, datadir, n_threads):
        '''--n-threads [N] --pipeline quick annotate --user-database [PEP.fa] [INPUT.fa]
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

    @pytest.mark.parametrize('n_threads', (1,4))
    def test_annotate_multiple_user_databases(self, tmpdir, datadir, n_threads):
        '''--pipeline quick annotate --user-database [PEP1.fa] --user-database [PEP2.fa] [INPUT.fa]
        '''

        with tmpdir.as_cwd():
            transcripts = datadir('pom.20.fa')
            pep = datadir('pep.fa')
            pep2 = datadir('odb_subset.fa')
            exp_gff3 = datadir('pom.20.udbs.dammit.gff3')
            exp_fasta = datadir('pom.20.udbs.dammit.fasta')

            args = ['run', '--n-threads', str(n_threads),
                    '--pipeline', 'quick', 'annotate',
                    '--user-database', pep,
                    '--user-database', pep2,
                    transcripts]
            status, out, err = run(*args)

            outdir = 'pom.20.dammit'
            gff3_fn = os.path.join(outdir, 'pom.20.dammit.gff3')
            fasta_fn = os.path.join(outdir, 'pom.20.dammit.fasta')

            assert status == 0
            assert compare_gff(gff3_fn, exp_gff3)
            assert open(fasta_fn).read() == open(exp_fasta).read()

    def test_annotate_basename(self, tmpdir, datadir):
        '''Test annotate --pipeline quick annotate --base-name [NAME] [INPUT.fa]
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

    def test_multiple_busco_groups(self, tmpdir, datadir):
        '''--pipeline quick --busco-group bacteria_odb10 --busco-group saccharomycetes_odb10
        '''

        with tmpdir.as_cwd():
            transcripts = datadir('pom.256.fa')
            exp_gff3 = datadir('pom.256.dammit.busco-multi.gff3')

            args = ['run', '--pipeline', 'quick',
                    '--busco-group', 'bacteria_odb10',
                    '--busco-group', 'saccharomycetes_odb10',
                    'annotate', transcripts]
            status, out, err = run(*args)
            assert status == 0

            gff3_fn = os.path.join('pom.256.dammit', 'pom.256.dammit.gff3')

            assert compare_gff(gff3_fn, exp_gff3)



    def test_annotate_outdir(self, tmpdir, datadir):
        '''
        Test output directory option
        '''

        with tmpdir.as_cwd():
            transcripts = datadir('pom.20.fa')

            outdir = 'test_out'
            args = ['run', '--pipeline', 'quick', 'annotate',
                    transcripts, '--output-dir', outdir]
            status, out, err = run(*args)
            assert os.path.isfile(os.path.join(outdir, 'pom.20.fasta'))

    # make sure DAMMIT_DB_DIR is set in your testing env
    # (export DAMMIT_DB_DIR=/path/to/databases)
    def test_annotate_dbdir_fail(self, tmpdir, datadir):
        '''Test annotation with a faulty database directory.
           dammit run --database-dir [DB_DIR] annotate [INPUT.fa]
        '''

        with tmpdir.as_cwd():
            transcripts = datadir('pom.20.fa')

            args = ['run', '--pipeline', 'quick', '--database-dir', '.', 'annotate', transcripts]
            status, out, err = run(*args, fail_ok=True)
            print(status, out, err)

            assert 'you probably need to install the dammit databases' in err
            assert status == 1


    def test_annotate_dbdir(self, tmpdir, datadir):
        '''Test that --database-dir works.
        '''

        with tmpdir.as_cwd():
            transcripts = datadir('pom.20.fa')
            database_dir = os.environ['DAMMIT_DB_DIR']
            exp_gff3 = datadir('pom.20.udb.dammit.gff3')
            exp_fasta = datadir('pom.20.udb.dammit.fasta')

            args = ['run', '--database-dir', database_dir, '--pipeline', 'quick', 'annotate',  transcripts]
            status, out, err = run(*args)

            outdir = 'pom.20.dammit'
            gff3_fn = os.path.join(outdir, 'pom.20.dammit.gff3')
            fasta_fn = os.path.join(outdir, 'pom.20.dammit.fasta')

            assert status == 0
            assert compare_gff(gff3_fn, exp_gff3)
            assert open(fasta_fn).read() == open(exp_fasta).read()


    def test_temp_dir(self, tmpdir, datadir):
        '''Test that --temp-dir works.
        '''

        with tmpdir.as_cwd():
            transcripts = datadir('pom.20.fa')
            dammit_temp_dir = "."
            args = ['run', '--temp-dir', dammit_temp_dir, '--pipeline', 'quick', 'annotate',  transcripts]
            status, out, err = run(*args)

            outdir = 'pom.20.dammit'

            assert status == 0
            tempd_contents = os.listdir(dammit_temp_dir)
            assert "pom.20.fa" in tempd_contents
            assert "pom.20.dammit" in tempd_contents


    def test_busco_group(self, tmpdir, datadir):
        '''Test that --busco-group works.
        '''

        with tmpdir.as_cwd():
            transcripts = datadir('pom.20.fa')
            dammit_temp_dir = "."
            args = ['run', '--busco-group', 'bacteria_odb10', '--pipeline', 'quick', 'annotate',  transcripts]
            status, out, err = run(*args)
            outdir = 'pom.20.dammit'

            assert status == 0
            assert os.path.isfile(os.path.join(outdir, "pom.20.busco/bacteria_odb10_outputs/run_bacteria_odb10/short_summary.txt"))

    def test_max_threads_per_task(self, tmpdir, datadir):
        '''Test that --max_threads_per_task works.
        '''

        with tmpdir.as_cwd():
            transcripts = datadir('pom.20.fa')
            args = ['run', '--max-threads-per-task', 1, '--pipeline', 'quick', 'annotate',  transcripts]
            status, out, err = run(*args)
            outdir = 'pom.20.dammit'

            print(status, out, err)

            assert status == 0
            assert "Threads (per-task):       1" in out

## do we not enable this anymore?
    #def test_config_file(self, tmpdir, datadir):
    #    '''Test that --config-file works.
    #    '''
#
 #       with tmpdir.as_cwd():
 #           transcripts = datadir('pom.20.fa')
 #           conf = datadir('test-conf.yml')
 #           args = ['run', '--config-file', conf, 'annotate',  transcripts]
 #           status, out, err = run(*args)
 #           outdir = 'pom.20.dammit'
#
#            print(status, out, err)
#
#            assert status == 0


    def test_busco_config_file(self, tmpdir, datadir):
        '''Test that --busco-config-file works.
        '''

        with tmpdir.as_cwd():
            transcripts = datadir('pom.20.fa')
            conf = datadir('test-busco-conf.ini')
            args = ['run', '--busco-config-file', conf, 'annotate',  transcripts]
            status, out, err = run(*args)
            outdir = 'pom.20.dammit'

            print(status, out, err)
            assert status == 0





