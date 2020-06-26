import os

from utils import run, runscript


def test_last(snakemake_rule, tmpdir):
    with tmpdir.as_cwd():
        snakemake_rule('last/last.rule', target='lastdb_transcript')
        
        assert os.path.isfile('test-transcript.fa.prj')
        assert os.path.isfile('test-transcript.log')
        assert 'lastdb: done!' in open('test-transcript.log').read()
