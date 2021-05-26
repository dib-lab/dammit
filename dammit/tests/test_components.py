# Tests for dammit CLI components related to FASTA munging:
# * dammit rename-fasta
# * dammit transcriptome-stats
# * dammit annotate-fasta

import filecmp
import json
import os
import yaml

from khmer import ReadParser
import pandas as pd

from .utils import run, runscript


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
    
    def test_basename_bad_chars(self, tmpdir, datadir):
        '''--basename with invalid characters fails'''
        with tmpdir.as_cwd():
            input_fa = datadir('test-transcript.fa')
            names_fn = 'names.csv'
            renamed_fn = 'renamed.fa'

            status, out, err = run('rename-fasta', '--basename', 'Test%$', input_fa, renamed_fn, names_fn)
            assert status == 1
            assert 'conform' in err

    def test_split_regex(self, tmpdir, datadir):
        '''--split-regex (?P<name>^[a-zA-Z0-9]+) should produce SPAC212'''

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
    
    def test_fail_on_repeats_fail(self, tmpdir, datadir):
        '''--fail-on-repeats should fail with exit code 1 on repeated header names'''
        with tmpdir.as_cwd():
            input_fa = datadir('repeated-names.fa')
            names_fn = 'names.csv'
            renamed_fn = 'renamed.fa'

            status, out, err = run('rename-fasta', '--fail-on-repeats',
                                   input_fa, renamed_fn, names_fn)
            assert status == 1
            assert 'ERROR: Repeated sequence names' in err

    def test_fail_on_repeats_pass(self, tmpdir, datadir):
        '''--fail-on-repeats should still pass in the absence of repeated header names'''
        with tmpdir.as_cwd():
            input_fa = datadir('pom.20.fa')
            names_fn = 'names.csv'
            renamed_fn = 'renamed.fa'

            status, out, err = run('rename-fasta', '--fail-on-repeats',
                                   input_fa, renamed_fn, names_fn)
            assert status == 0
            
            names = pd.read_csv(names_fn)

            # OOPS pom.20.fa actually has 21 sequences loooooool
            assert len(names) == 21


class TestTranscriptomeStats:
    '''dammit transcriptome-stats'''

    def test_defaults(self, tmpdir, datadir):
        '''defaults (k=25) produce correct numbers'''
        with tmpdir.as_cwd():
            transcript = datadir('test-transcript.fa')
            output_fn = str(tmpdir.join('test'))
            lens_fn = str(tmpdir.join('lens'))

            status, out, err = run('transcriptome-stats', transcript, output_fn, lens_fn)

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
            lens_fn = str(tmpdir.join('lens'))

            status, out, err = run('transcriptome-stats', '-K', '27', transcript, output_fn, lens_fn)

            with open(output_fn) as fp:
                results = json.load(fp)

            assert results['k_mers'] == 5636
            assert results['k_mers_unique'] == 5627

    def test_contains_non_ACGTN(self, tmpdir, datadir):
        '''non-ACGTN characters throw a RunTimeError'''

        with tmpdir.as_cwd():
            transcript = datadir('non-actg-transcripts.fa')
            output_fn = str(tmpdir.join('test'))
            lens_fn = str(tmpdir.join('lens'))

            status, out, err = run('transcriptome-stats', transcript, output_fn, lens_fn)
            
            assert status == -1
            assert 'Offending' in err

    def test_contains_N(self, tmpdir, datadir):
        '''counts N's in transcripts properly'''
        with tmpdir.as_cwd():
            transcript = datadir('test-transcript-N.fa')
            output_fn = str(tmpdir.join('test'))
            lens_fn = str(tmpdir.join('lens'))

            status, out, err = run('transcriptome-stats', transcript, output_fn, lens_fn)
            
            with open(output_fn) as fp:
                results = json.load(fp)

            assert 'n_ambiguous' in results
            assert results['n_ambiguous'] == 10

            assert round(results['GCperc'], 4) == 0.4152


class TestAnnotateFasta:
    '''dammit annotate-fasta'''

    def test_defaults(self, tmpdir, datadir):
        '''annotates (renamed) pom.single.fa from provided GFF3 file'''
        with tmpdir.as_cwd():
            base_fa = datadir('pom.single.fa')
            renamed_fa = 'pom.renamed.fa'
            gff3 = datadir('pom.single.fa.dammit.gff3')
            output_fa = 'pom.annotated.fa'

            run('rename-fasta', base_fa, renamed_fa, 'names.csv')
            run('annotate-fasta', renamed_fa, gff3, output_fa)

            sequence_header = list(ReadParser(output_fa))[0].name
            assert sequence_header == 'Transcript_0 len=5662 CDS=0-5661 exon=0-5662 gene=0-5662 mRNA=0-5662 hmm_matches=DEAD:3603-4089,Helicase_C:4206-4548,Helicase_C:5217-5274 three_prime_UTR=5661-5662 homologies=TLH2_SCHPO:0-5661'
    
    def test_name_map(self, tmpdir, datadir):
        '''annotates pom.single.fa from provided GFF3 file and maps back to original names.'''
        with tmpdir.as_cwd():
            base_fa = datadir('pom.single.fa')
            renamed_fa = 'pom.renamed.fa'
            gff3 = datadir('pom.single.fa.dammit.gff3')
            output_fa = 'pom.annotated.fa'
            expected_fa = datadir('pom.single.annotated.fa')

            run('rename-fasta', base_fa, renamed_fa, 'names.csv')
            run('annotate-fasta', renamed_fa, gff3, output_fa, '--name-map', 'names.csv')

            sequence_header = list(ReadParser(output_fa))[0].name
            expected_header = list(ReadParser(expected_fa))[0].name
            assert sequence_header == expected_header


class TestBestHits:
    '''dammit best-hits'''

    def test_defaults(self, tmpdir, datadir):
        '''extracts the top scoring hit for pom.single.fa against OrthoDB'''
        with tmpdir.as_cwd():
            input_maf = datadir('pom.single.fa.x.OrthoDB.maf')
            output_csv = 'best.csv'

            run('best-hits', input_maf, output_csv)

            hits_df = pd.read_csv(output_csv)
            assert hits_df.s_name[0] == 'TLH2_SCHPO'
            assert hits_df.E[0] == 0
            assert hits_df.score[0] == 10676.0
            assert hits_df.s_len[0] == 1919


class TestMafToGFF3:
    '''dammit maf-to-gff3'''

    def check_pom_gff3(self, filename):
        raw_gff3 = open(filename).readlines()
        assert len(raw_gff3) == 47

        row = raw_gff3[1].split('\t')
        assert row[0] == 'Transcript_0'
        assert row[3] == '1'
        assert row[4] == '5661'
        assert row[6] == '+'
        assert 'Target=TLH2_SCHPO 1 1887' in row[-1]
        assert ';database=OrthoDB' in row[-1]

        return row

    def test_maf_input(self, tmpdir, datadir):
        '''converts raw MAF input to corresponding GFF3'''
        with tmpdir.as_cwd():
            input_maf = datadir('pom.single.fa.x.OrthoDB.maf')
            output_gff3 = 'maf.gff3'

            run('maf-to-gff3', '--dbxref', 'OrthoDB', input_maf, output_gff3)

            row = self.check_pom_gff3(output_gff3)
            assert 'homology:58381196ee8824e0ae85ca41d72cabb7779270b5' in row[-1]
            assert row[2] == 'translated_nucleotide_match'

    def test_csv_input(self, tmpdir, datadir):
        '''converts shmlast-style MAF CSV input to GFF3'''
        with tmpdir.as_cwd():
            input_csv = datadir('pom.single.fa.x.OrthoDB.maf.csv')
            output_gff3 = 'maf.gff3'

            run('maf-to-gff3', '--dbxref', 'OrthoDB', input_csv, output_gff3)

            row = self.check_pom_gff3(output_gff3)
            assert 'homology:58381196ee8824e0ae85ca41d72cabb7779270b5' in row[-1]
            assert row[2] == 'translated_nucleotide_match'


    def test_shmlast_input(self, tmpdir, datadir):
        '''works the same as shmlast-to-gff3'''
        with tmpdir.as_cwd():
            input_csv = datadir('pom.single.fa.x.OrthoDB.maf.csv')
            output_gff3 = 'maf.gff3'

            run('shmlast-to-gff3', '--dbxref', 'OrthoDB', input_csv, output_gff3)

            row = self.check_pom_gff3(output_gff3)
            assert row[1] == 'shmlast.LAST'
            assert row[2] == 'conditional_reciprocal_best_LAST'


class TestHMMERToGFF3:
    '''dammit hmmscan-to-gff3'''

    def test_defaults(self, tmpdir, datadir):
        '''converts hmmscan tbl to proper GFF3'''
        with tmpdir.as_cwd():
            input_tbl = datadir('test-protein-x-pfam-a.tbl')
            output_gff3 = 'protein.gff3'
            expected_gff3 = datadir('test-protein-x-pfam-a.gff3')

            run('hmmscan-to-gff3', '--dbxref', 'Pfam', input_tbl, output_gff3)

            assert filecmp.cmp(output_gff3, expected_gff3, shallow=False)


class TestInfernalToGFF3:
    '''dammit cmscan-to-gff3'''

    def test_defaults(self, tmpdir, datadir):
        '''converts cmscan tbl to proper GFF3'''
        with tmpdir.as_cwd():
            input_tbl = datadir('rnaseP.tbl')
            output_gff3 = 'rnaseP.out.gff3'
            expected_gff3 = datadir('rnaseP.tbl.gff3')

            run('cmscan-to-gff3', '--dbxref', 'Rfam', input_tbl, output_gff3)

            assert filecmp.cmp(output_gff3, expected_gff3, shallow=False)


class TestMergeGFF3:
    '''dammit merge-gff3'''

    def test_defaults(self, tmpdir, datadir):
        '''merged gff3 files are sane'''
        with tmpdir.as_cwd():
            input_gff3_1 = datadir('test-protein-x-pfam-a.gff3')
            input_gff3_2 = datadir('rnaseP.tbl.gff3')
            expected = datadir('merged.gff3')
            actual = 'actual.gff3'

            run('merge-gff3', input_gff3_1, input_gff3_2, actual)

            assert filecmp.cmp(expected, actual, shallow=False)

    def test_ordering(self, tmpdir, datadir):
        '''merged gff3 files are agnostic to input order'''
        with tmpdir.as_cwd():
            input_gff3_1 = datadir('test-protein-x-pfam-a.gff3')
            input_gff3_2 = datadir('rnaseP.tbl.gff3')
            expected = datadir('merged.gff3')
            actual = 'actual.gff3'

            run('merge-gff3', input_gff3_2, input_gff3_1, actual)

            assert filecmp.cmp(expected, actual, shallow=False)


class TestConfigInfo:
    '''dammit config'''

    def test_show_default_core(self, tmpdir):
        '''show-defaults core produces reasonable yaml'''

        status, out, err = run('config', 'show-default', 'core')
        assert status == 0

        data = yaml.safe_load(out)

        assert data['basename'] == 'Transcript'
        assert data['output_dir'] == None
        assert data['dammit_dir_suffix'] == '.dammit'
        assert data['input_transcriptome'] == None

    def test_show_default_databases(self, tmpdir):
        '''show-defaults databases produces reasonable yaml'''

        status, out, err = run('config', 'show-default', 'databases')
        assert status == 0

        data = yaml.safe_load(out)

        assert set(['OrthoDB', 'Pfam-A', 'Rfam', 'busco', 'busco', 'nr']) < set(data.keys())

    def test_show_default_pipelines(self, tmpdir):
        '''show-defaults pipelines produces reasonable yaml'''

        status, out, err = run('config', 'show-default', 'pipelines')
        assert status == 0

        data = yaml.safe_load(out)

        assert set(['default', 'full', 'nr', 'quick']) <= set(data['pipelines'].keys())
    
    def test_show_directories(self):
        '''show-directories contains the right keys'''
        status, out, err = run('config', 'show-directories')
        assert status == 0

        for item in ['Databases dir:', 'Temp dir:', 'Conda env dir:']:
            assert item in out
    
    def test_busco_groups(self):
        '''busco-groups contains enough entries'''
        status, out, err = run('config', 'busco-groups')
        assert status == 0
        
        groups = out.split(' ')
        assert 'bacteria_odb10' in groups
        assert 'solanales_odb10' in groups
        