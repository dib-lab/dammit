# Tests for dammit CLI components related to FASTA munging:
# * dammit rename-fasta
# * dammit transcriptome-stats
# * dammit annotate-fasta

import json
import os

from khmer import ReadParser
import pandas as pd

from utils import run, runscript

class TestRenameFasta:
    '''dammit rename-fasta'''

    def test_defaults(self, tmpdir, datadir):
        '''defaults should produce Transcript_0'''

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


    def test_basename(self, tmpdir, datadir):
        '''--basename Test should produce Test_0'''
        with tmpdir.as_cwd():
            input_fa = datadir('test-transcript.fa')
            names_fn = 'names.csv'
            renamed_fn = 'renamed.fa'

            run('rename-fasta', '--basename', 'Test', input_fa, renamed_fn, names_fn)

            names_df = pd.read_csv(names_fn)
            assert names_df['renamed'][0] == 'Test_0'

            renamed_records = list(ReadParser(renamed_fn))
            assert renamed_records[0].name == 'Test_0'


    def test_split_regex(self, tmpdir, datadir):
        '''--split-refex (?P<name>^[a-zA-Z0-9]+) should produce SPAC212'''

        with tmpdir.as_cwd():
            input_fa = datadir('test-transcript.fa')
            names_fn = 'names.csv'
            renamed_fn = 'renamed.fa'


            run('rename-fasta', '--split-regex', r'(?P<name>^[a-zA-Z0-9]+)', input_fa, renamed_fn, names_fn)

            names_df = pd.read_csv(names_fn)
            assert names_df['renamed'][0] == 'SPAC212'

            renamed_records = list(ReadParser(renamed_fn))
            assert renamed_records[0].name == 'SPAC212'


    def test_conflicting(self, tmpdir, datadir):
        '''--split-refex (?P<name>^[a-zA-Z0-9]+) should override --basename Test'''
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


class TestTranscriptomeStats:
    '''dammit transcriptome-stats'''

    def test_defaults(self, tmpdir, datadir):
        '''defaults (k=25) produce correct numbers'''
        with tmpdir.as_cwd():
            transcript = datadir('test-transcript.fa')
            output_fn = str(tmpdir.join('test'))

            status, out, err = run('transcriptome-stats', transcript, output_fn)

            with open(output_fn) as fp:
                results = json.load(fp)

            assert 'n_ambiguous' in results
            assert results['n_ambiguous'] == 0

            assert 'N' in results
            assert results['N'] == 1

            assert results['k_mers'] == 5638
            assert results['k_mers_unique'] == 5653

            assert round(results['GCperc'], 4) == 0.4149
    
    def test_ksize(self, tmpdir, datadir):
        '''--K 27 produces correct k-mer counts'''
        with tmpdir.as_cwd():
            transcript = datadir('test-transcript.fa')
            output_fn = str(tmpdir.join('test'))

            status, out, err = run('transcriptome-stats', '-K', '27', transcript, output_fn)

            with open(output_fn) as fp:
                results = json.load(fp)

            assert results['k_mers'] == 5636
            assert results['k_mers_unique'] == 5627

    def test_contains_non_ACGTN(self, tmpdir, datadir):
        '''non-ACGTN characters throw a RunTimeError'''

        with tmpdir.as_cwd():
            transcript = datadir('non-actg-transcripts.fa')
            output_fn = str(tmpdir.join('test'))

            status, out, err = run('transcriptome-stats', transcript, output_fn)
            
            assert status == -1
            assert 'Offending' in err

    def test_contains_N(self, tmpdir, datadir):
        '''counts N's in transcripts properly'''
        with tmpdir.as_cwd():
            transcript = datadir('test-transcript-N.fa')
            output_fn = str(tmpdir.join('test'))

            status, out, err = run('transcriptome-stats', transcript, output_fn)
            
            with open(output_fn) as fp:
                results = json.load(fp)

            assert 'n_ambiguous' in results
            assert results['n_ambiguous'] == 10

            assert round(results['GCperc'], 4) == 0.4152