import os

from utils import run, runscript


def test_lastdb_transcript(snakemake_rule, tmpdir):
    with tmpdir.as_cwd():
        snakemake_rule('last/last.rule', target='lastdb_transcript')

        assert os.path.isfile('test-transcript.fa.prj')
        assert os.path.isfile('test-transcript.log')
        assert 'lastdb: done!' in open('test-transcript.log').read()


def test_lastdb_protein(snakemake_rule, tmpdir):
    with tmpdir.as_cwd():
        snakemake_rule('last/last.rule', target='lastdb_protein')

        assert os.path.isfile('test-protein.fa.prj')
        assert os.path.isfile('test-protein.log')
        assert 'lastdb: done!' in open('test-protein.log').read()


def test_lastal_nucl_x_nucl(snakemake_rule, tmpdir):
    with tmpdir.as_cwd():
        snakemake_rule('last/last.rule', target='lastal_nucl_x_nucl')

        assert os.path.isfile('test-transcript.maf')
        assert os.path.isfile('lastal_nucl_x_nucl.log')
        assert 'test-transcript.fa' in open('test-transcript.maf').read()


def test_lastal_nucl_x_prot(snakemake_rule, tmpdir):
    with tmpdir.as_cwd():
        snakemake_rule('last/last.rule', target='lastal_nucl_x_prot')
        assert os.path.isfile('test-tr-x-prot.maf')
        assert os.path.isfile('lastal_nucl_x_prot.log')
        assert 'test-protein.fa' in open('test-tr-x-prot.maf').read()



