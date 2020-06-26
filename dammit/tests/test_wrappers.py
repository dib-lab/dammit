import os
import pytest
from utils import run, runscript


def test_lastdb_transcript_dryrun(snakemake_rule, tmpdir):
    with tmpdir.as_cwd():
        status, out, err = snakemake_rule('last/last.rule', target='lastdb_transcript', extra_args =["-n"])

        assert status == 0


@pytest.mark.long
def test_lastdb_transcript(snakemake_rule, tmpdir):
    with tmpdir.as_cwd():
        snakemake_rule('last/last.rule', target='lastdb_transcript')

        assert os.path.isfile('test-transcript.fa.prj')
        assert os.path.isfile('test-transcript.log')
        assert 'lastdb: done!' in open('test-transcript.log').read()


def test_lastdb_protein_dryrun(snakemake_rule, tmpdir):
    with tmpdir.as_cwd():
        status, out, err = snakemake_rule('last/last.rule', target='lastdb_protein', extra_args =["-n"])

        assert status == 0


@pytest.mark.long
def test_lastdb_protein(snakemake_rule, tmpdir):
    with tmpdir.as_cwd():
        snakemake_rule('last/last.rule', target='lastdb_protein')

        assert os.path.isfile('test-protein.fa.prj')
        assert os.path.isfile('test-protein.log')
        assert 'lastdb: done!' in open('test-protein.log').read()


def test_lastal_nucl_x_nucl_dryrun(snakemake_rule, tmpdir):
    with tmpdir.as_cwd():
        status, out, err = snakemake_rule('last/last.rule', target='lastal_nucl_x_nucl', extra_args =["-n"])

        assert status == 0


@pytest.mark.long
def test_lastal_nucl_x_nucl(snakemake_rule, tmpdir):
    with tmpdir.as_cwd():
        snakemake_rule('last/last.rule', target='lastal_nucl_x_nucl')

        assert os.path.isfile('test-transcript.maf')
        assert os.path.isfile('lastal_nucl_x_nucl.log')
        assert 'test-transcript.fa' in open('test-transcript.maf').read()


def test_lastal_nucl_x_prot_dryrun(snakemake_rule, tmpdir):
    with tmpdir.as_cwd():
        status, out, err = snakemake_rule('last/last.rule', target='lastal_nucl_x_prot', extra_args =["-n"])

        assert status == 0


@pytest.mark.long
def test_lastal_nucl_x_prot(snakemake_rule, tmpdir):
    with tmpdir.as_cwd():
        snakemake_rule('last/last.rule', target='lastal_nucl_x_prot')

        assert os.path.isfile('test-tr-x-prot.maf')
        assert os.path.isfile('lastal_nucl_x_prot.log')
        assert 'test-protein.fa' in open('test-tr-x-prot.maf').read()


def test_transdecoder_longorfs_dryrun(snakemake_rule, tmpdir):
    with tmpdir.as_cwd():
        status, out, err = snakemake_rule('transdecoder/transdecoder.rule', target='transdecoder_longorfs', extra_args =["-n"])

        assert status == 0


@pytest.mark.long
def test_transdecoder_longorfs(snakemake_rule, tmpdir):
    with tmpdir.as_cwd():
        snakemake_rule('transdecoder/transdecoder.rule', target='transdecoder_longorfs')

        assert os.path.isfile('test.fa.transdecoder_dir/longest_orfs.pep')
        assert os.path.isfile('test-longorfs.log')
        assert 'Done preparing long ORFs' in open('test-longorfs.log').read()


def test_transdecoder_predict_dryrun(snakemake_rule, tmpdir):
    with tmpdir.as_cwd():
        status, out, err = snakemake_rule('transdecoder/transdecoder.rule', target='transdecoder_predict', extra_args =["-n"])

        assert status == 0

@pytest.mark.long
def test_transdecoder_predict(snakemake_rule, tmpdir):
    with tmpdir.as_cwd():
        snakemake_rule('transdecoder/transdecoder.rule', target='transdecoder_predict')

        assert os.path.isfile('test.fa.transdecoder.gff3')
        assert os.path.isfile('test-predict.log')
        assert 'transdecoder is finished' in open('test-predict.log').read()

