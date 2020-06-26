# Tests for dammit CLI components related to FASTA munging:
# * dammit rename-fasta
# * dammit transcriptome-stats
# * dammit annotate-fasta


import os

from khmer import ReadParser
import pandas as pd

from utils import run, runscript


def test_rename_fasta_defaults(tmpdir, datadir):
    with tmpdir.as_cwd():
        input_fa = datadir('test-transcript.fa')
        names_fn = 'names.csv'
        renamed_fn = 'renamed.fa'

        run('rename-fasta', input_fa, renamed_fn, names_fn)

        assert os.path.isfile(names_fn)
        assert os.path.isfile(renamed_fn)

        names_df = pd.read_csv(names_fn)
        assert names_df['original'][0] == 'SPAC212_RecQ_type_DNA_helicase_TRANSCRIPT'
        assert names_df['renamed'][0] == 'Transcript_0'

        original_records = list(ReadParser(input_fa))
        renamed_records = list(ReadParser(renamed_fn))

        assert original_records[0].sequence == renamed_records[0].sequence
        assert renamed_records[0].name == 'Transcript_0'


def test_rename_fasta_basename(tmpdir, datadir):
    with tmpdir.as_cwd():
        input_fa = datadir('test-transcript.fa')
        names_fn = 'names.csv'
        renamed_fn = 'renamed.fa'

        run('rename-fasta', '--basename', 'Test', input_fa, renamed_fn, names_fn)

        names_df = pd.read_csv(names_fn)
        assert names_df['renamed'][0] == 'Test_0'

        renamed_records = list(ReadParser(renamed_fn))
        assert renamed_records[0].name == 'Test_0'


def test_rename_fasta_regex(tmpdir, datadir):
    with tmpdir.as_cwd():
        input_fa = datadir('test-transcript.fa')
        names_fn = 'names.csv'
        renamed_fn = 'renamed.fa'


        run('rename-fasta', '--split-regex', r'(?P<name>^[a-zA-Z0-9]+)', input_fa, renamed_fn, names_fn)

        names_df = pd.read_csv(names_fn)
        assert names_df['renamed'][0] == 'SPAC212'

        renamed_records = list(ReadParser(renamed_fn))
        assert renamed_records[0].name == 'SPAC212'


def test_rename_fasta_conflicting(tmpdir, datadir):
    with tmpdir.as_cwd():
        input_fa = datadir('test-transcript.fa')
        names_fn = 'names.csv'
        renamed_fn = 'renamed.fa'


        status, out, err = run('rename-fasta', '--basename', 'Test',
                               '--split-regex', r'(?P<name>^[a-zA-Z0-9]+)', 
                               input_fa, renamed_fn, names_fn)
        assert 'NOTE:' in err

        names_df = pd.read_csv(names_fn)
        assert names_df['renamed'][0] == 'SPAC212'

        renamed_records = list(ReadParser(renamed_fn))
        assert renamed_records[0].name == 'SPAC212'