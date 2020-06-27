import os

from utils import run, runscript


def test_last(snakemake_rule, tmpdir, datadir):
    with tmpdir.as_cwd():
        input_fn = datadir('test-transcript.fa')
        status, out, err = snakemake_rule('last/last.rule',
                                          target='lastdb_transcript',
                                          config={'DATA_DIR': str(tmpdir)})
        
        assert os.path.isfile('test-transcript.fa.prj')
        assert os.path.isfile('test-transcript.log')
        assert 'lastdb: done!' in open('test-transcript.log').read()