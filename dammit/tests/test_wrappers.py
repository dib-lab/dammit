import os
import pytest
from .utils import run, runscript


def test_lastdb_transcript_dryrun(snakemake_rule, tmpdir, datadir):
    with tmpdir.as_cwd():
        datadir('test-transcript.fa')
        status, out, err = snakemake_rule('last/last.rule',
                                          target='lastdb_transcript',
                                          config={'DATA_DIR': str(tmpdir)},
                                          extra_args =["-n"])

        assert status == 0


@pytest.mark.long
def test_lastdb_transcript(snakemake_rule, tmpdir, datadir):
    with tmpdir.as_cwd():
        datadir('test-transcript.fa')
        status, out, err = snakemake_rule('last/last.rule',
                                          target='lastdb_transcript',
                                          config={'DATA_DIR': str(tmpdir)})

        assert os.path.isfile('test-transcript.fa.prj')
        assert os.path.isfile('test-transcript.log')
        assert 'lastdb: done!' in open('test-transcript.log').read()


def test_lastdb_protein_dryrun(snakemake_rule, tmpdir, datadir):
    with tmpdir.as_cwd():
        datadir('test-protein.fa')
        status, out, err = snakemake_rule('last/last.rule',
                                           target='lastdb_protein',
                                           config={'DATA_DIR': str(tmpdir)},
                                           extra_args =["-n"])

        assert status == 0


@pytest.mark.long
def test_lastdb_protein(snakemake_rule, tmpdir, datadir):
    with tmpdir.as_cwd():
        datadir('test-protein.fa')
        status, out, err = snakemake_rule('last/last.rule',
                                           target='lastdb_protein',
                                           config={'DATA_DIR': str(tmpdir)})

        assert os.path.isfile('test-protein.fa.prj')
        assert os.path.isfile('test-protein.log')
        assert 'lastdb: done!' in open('test-protein.log').read()


def test_lastal_nucl_x_nucl_dryrun(snakemake_rule, tmpdir, datadir):
    with tmpdir.as_cwd():
        datadir('test-transcript.fa')
        status, out, err = snakemake_rule('last/last.rule', \
                                           target='lastal_nucl_x_nucl',
                                           config={'DATA_DIR': str(tmpdir)},
                                           extra_args =["-n"])

        assert status == 0


@pytest.mark.long
def test_lastal_nucl_x_nucl(snakemake_rule, tmpdir, datadir):
    with tmpdir.as_cwd():
        datadir('test-transcript.fa')
        status, out, err = snakemake_rule('last/last.rule',
                                           target='lastal_nucl_x_nucl',
                                           config={'DATA_DIR': str(tmpdir)})


        assert os.path.isfile('test-transcript.maf')
        assert os.path.isfile('lastal_nucl_x_nucl.log')
        assert 'test-transcript.fa' in open('test-transcript.maf').read()


def test_lastal_nucl_x_prot_dryrun(snakemake_rule, tmpdir, datadir):
    with tmpdir.as_cwd():
        datadir('test-transcript.fa')
        datadir('test-protein.fa')
        status, out, err = snakemake_rule('last/last.rule',
                                           target='lastal_nucl_x_prot',
                                           config={'DATA_DIR': str(tmpdir)},
                                           extra_args =["-n"])

        assert status == 0


@pytest.mark.long
def test_lastal_nucl_x_prot(snakemake_rule, tmpdir, datadir):
    with tmpdir.as_cwd():
        datadir('test-transcript.fa')
        datadir('test-protein.fa')
        status, out, err = snakemake_rule('last/last.rule',
                                           target='lastal_nucl_x_prot',
                                           config={'DATA_DIR': str(tmpdir)})

        assert os.path.isfile('test-tr-x-prot.maf')
        assert os.path.isfile('lastal_nucl_x_prot.log')
        assert 'test-protein.fa' in open('test-tr-x-prot.maf').read()


def test_transdecoder_longorfs_dryrun(snakemake_rule, tmpdir, datadir):
    with tmpdir.as_cwd():
        datadir('test.fa.gz')
        status, out, err = snakemake_rule('transdecoder/transdecoder.rule',
                                           target='transdecoder_longorfs',
                                           config={'DATA_DIR': str(tmpdir)},
                                           extra_args =["-n"])

        assert status == 0


@pytest.mark.long
def test_transdecoder_longorfs(snakemake_rule, tmpdir, datadir):
    with tmpdir.as_cwd():
        datadir('test.fa.gz')
        status, out, err = snakemake_rule('transdecoder/transdecoder.rule',
                                           target='transdecoder_longorfs',
                                           config={'DATA_DIR': str(tmpdir)})


        assert os.path.isfile('test.fa.transdecoder_dir/longest_orfs.pep')
        assert os.path.isfile('test-longorfs.log')
        assert 'Done preparing long ORFs' in open('test-longorfs.log').read()


def test_transdecoder_predict_dryrun(snakemake_rule, tmpdir, datadir):
    with tmpdir.as_cwd():
        datadir('test.fa.gz')
        status, out, err = snakemake_rule('transdecoder/transdecoder.rule',
                                           target='transdecoder_predict',
                                           config={'DATA_DIR': str(tmpdir)},
                                           extra_args =["-n"])

        assert status == 0

@pytest.mark.long
def test_transdecoder_predict(snakemake_rule, tmpdir, datadir):
    with tmpdir.as_cwd():
        datadir('test.fa.gz')
        status, out, err = snakemake_rule('transdecoder/transdecoder.rule',
                                           target='transdecoder_predict',
                                           config={'DATA_DIR': str(tmpdir)})

        assert os.path.isfile('test.fa.transdecoder.gff3')
        assert os.path.isfile('test-predict.log')
        assert 'transdecoder is finished' in open('test-predict.log').read()


def test_infernal_cmpress_dryrun(snakemake_rule, tmpdir, datadir):
    with tmpdir.as_cwd():
        datadir('test-covariance-model.cm')
        status, out, err = snakemake_rule('infernal/infernal.rule',
                                           target='infernal_cmpress',
                                           config={'DATA_DIR': str(tmpdir)},
                                           extra_args=["-n"])

        assert status == 0


@pytest.mark.long
def test_infernal_cmpress(snakemake_rule, tmpdir, datadir):
    with tmpdir.as_cwd():
        datadir('test-covariance-model.cm')
        status, out, err = snakemake_rule('infernal/infernal.rule',
                                           target='infernal_cmpress',
                                           config={'DATA_DIR': str(tmpdir)})

        print(out)
        assert os.path.isfile('test-covariance-model.cm.i1i')
        assert os.path.isfile('test-cmpress.log')
        assert 'Working...    done.' in open('test-cmpress.log').read()


def test_infernal_cmscan_dryrun(snakemake_rule, tmpdir, datadir):
    with tmpdir.as_cwd():
        datadir('test-covariance-model.cm')
        datadir('test-transcript.fa')
        status, out, err = snakemake_rule('infernal/infernal.rule',
                                           target='infernal_cmscan',
                                           config={'DATA_DIR': str(tmpdir)},
                                           extra_args=["-n"])

        assert status == 0


@pytest.mark.long
def test_infernal_cmscan(snakemake_rule, tmpdir, datadir):
    with tmpdir.as_cwd():
        datadir('test-covariance-model.cm')
        datadir('test-transcript.fa')
        status, out, err = snakemake_rule('infernal/infernal.rule',
                                           target='infernal_cmscan',
                                           config={'DATA_DIR': str(tmpdir)})

        print(out)
        assert os.path.isfile('tr-infernal-tblout.txt')
        assert os.path.isfile('test-cmscan.log')
        assert '[ok]' in open('test-cmscan.log').read()

def test_busco_dryrun(snakemake_rule, tmpdir, datadir):
    with tmpdir.as_cwd():
        datadir('target.fa')
        datadir('busco-config.ini')
        status, out, err = snakemake_rule('busco/busco.rule',
                                           target='run_busco',
                                           config={'DATA_DIR': str(tmpdir)},
                                           extra_args=["-n"])

        assert status == 0


@pytest.mark.long
def test_busco(snakemake_rule, tmpdir, datadir):
    with tmpdir.as_cwd():
        datadir('target.fa')
        datadir('busco-config.ini')
        status, out, err = snakemake_rule('busco/busco.rule',
                                           target='run_busco',
                                           config={'DATA_DIR': str(tmpdir)})

        print(out)
        #assert os.path.isfile('txome_busco/full_table_txome_busco.tsv')
        assert os.path.isfile('test-busco.log')
        print(open('test-busco.log').read())
        assert 'Results from dataset metazoa_odb10' in open('test-busco.log').read()
        assert "C:0.1%[S:0.1%,D:0.0%],F:0.0%,M:99.9%,n:954" in open('test-busco.log').read()
        assert 'BUSCO analysis done.' in open('test-busco.log').read()


# TODO: explicitly test busco configuration?

def test_hmmscan_dryrun(snakemake_rule, tmpdir, datadir):
    with tmpdir.as_cwd():
        datadir('test-protein-multi.fa')
        #datadir('test-protein-multi.fa')
        datadir('test-profile.sto')
        status, out, err = snakemake_rule('hmmer/hmmer.rule',
                                           target='hmmscan_profile',
                                           config={'DATA_DIR': str(tmpdir)},
                                           extra_args=["-n"])

        assert status == 0


@pytest.mark.parametrize('n_threads', (1,4))
def test_hmmscan(snakemake_rule, tmpdir, datadir, n_threads):
    with tmpdir.as_cwd():
        datadir('test-protein-multi.fa')
        datadir('test-profile.hmm')
        status, out, err = snakemake_rule('hmmer/hmmer.rule',
                                           target='hmmscan_profile',
                                           config={'DATA_DIR': str(tmpdir)},
                                           n_threads=n_threads)

        assert status == 0
        print('STDOUT:', out)
        print('STDERR:', err)
        # check hmmpress
        assert "Pressed and indexed 1 HMMs (1 names and 1 accessions)." in open('hmmpress.log').read()
        assert "Profiles (remainder) pressed into: test-profile.hmm.h3p" in open('hmmpress.log').read()
        # check hmmscan
        assert os.path.isfile('hmmscan.log')
        ### not getting any hits in here; prior test on same data (in test_hmmer.py found 2 hits)
        print(open('hmmscan/test-prot-tbl.txt').read())
        print(open('hmmscan/test-prot-out.txt').read())


def test_hmmpress_dryrun(snakemake_rule, tmpdir, datadir):
    with tmpdir.as_cwd():
        datadir('test-profile.hmm')
        status, out, err = snakemake_rule('hmmer/hmmer.rule',
                                           target='hmmpress_profile',
                                           config={'DATA_DIR': str(tmpdir)},
                                           extra_args=["-n"])

        assert status == 0

def test_hmmpress(snakemake_rule, tmpdir, datadir):
    with tmpdir.as_cwd():
        datadir('test-profile.hmm')
        status, out, err = snakemake_rule('hmmer/hmmer.rule',
                                           target='hmmpress_profile',
                                           config={'DATA_DIR': str(tmpdir)})

        assert os.path.isfile('hmmpress.log')
        # check hmmpress
        assert "Pressed and indexed 1 HMMs (1 names and 1 accessions)." in open('hmmpress.log').read()
        assert "Profiles (remainder) pressed into: test-profile.hmm.h3p" in open('hmmpress.log').read()


def test_hmmbuild_dryrun(snakemake_rule, tmpdir, datadir):
    with tmpdir.as_cwd():
        datadir('test-profile.sto')
        status, out, err = snakemake_rule('hmmer/hmmer.rule',
                                           target='hmmbuild_profile',
                                           config={'DATA_DIR': str(tmpdir)},
                                           extra_args=["-n"])

        assert status == 0

def test_hmmbuild(snakemake_rule, tmpdir, datadir):
    with tmpdir.as_cwd():
        datadir('test-profile.sto')
        status, out, err = snakemake_rule('hmmer/hmmer.rule',
                                           target='hmmbuild_profile',
                                           config={'DATA_DIR': str(tmpdir)})

        assert os.path.isfile('test-profile-hmmbuild.log')
        #check hmmbuild
        assert all (("1     Ribosomal_L22        11354   603   325" in open("test-profile-hmmbuild.log").read(),
                   "17.71  0.590 Ribosomal protein L22p/L17e" in open("test-profile-hmmbuild.log").read()))


def test_shmlast_crbl_dryrun(snakemake_rule, tmpdir, datadir):
    with tmpdir.as_cwd():
        datadir('test-transcript.fa')
        datadir('pep.reduced.fa')
        status, out, err = snakemake_rule('shmlast/shmlast.rule',
                                           target='shmlast_crbl',
                                           config={'DATA_DIR': str(tmpdir)},
                                           extra_args=["-n"])

        assert status == 0


def test_shmlast_crbl(snakemake_rule, tmpdir, datadir):
    with tmpdir.as_cwd():
        datadir('test-transcript.fa')
        datadir('pep.reduced.fa')
        status, out, err = snakemake_rule('shmlast/shmlast.rule',
                                           target='shmlast_crbl',
                                           config={'DATA_DIR': str(tmpdir)})

        assert status == 0
        assert os.path.isfile('test-shmlast-crbl.log')
        assert "--- Begin Task Execution ---" in open("test-shmlast-crbl.log").read()
        assert all(("2.7e-88,2.5e-76,87.56863623584101,0,310.9231784374275,158,0,1888" in open("transcripts.x.pep.reduced.shmlast_crbl.csv").read(),
                    "SPAC212_RecQ_type_DNA_helicase_TRANSCRIPT" in open("transcripts.x.pep.reduced.shmlast_crbl.csv").read()))


def test_download_and_gunzip_dryrun(snakemake_rule, tmpdir, datadir):
    with tmpdir.as_cwd():
        datadir('test-transcript.fa')
        status, out, err = snakemake_rule('download/download.rule',
                                           target='download_and_gunzip_url',
                                           config={'DATA_DIR': str(tmpdir)},
                                           extra_args=["-n"])

        assert status == 0

def test_download_and_gunzip_url(snakemake_rule, tmpdir, datadir):
    with tmpdir.as_cwd():
        datadir('test-transcript.fa')
        status, out, err = snakemake_rule('download/download.rule',
                                           target='download_and_gunzip_url')
                                           #config={'DATA_DIR': str(tmpdir)})

        assert status == 0

def test_download_and_gunzip_url_md5(snakemake_rule, tmpdir, datadir):
    with tmpdir.as_cwd():
        datadir('test-transcript.fa')
        status, out, err = snakemake_rule('download/download.rule',
                                           target='download_and_gunzip_urL_md5')
                                           #config={'DATA_DIR': str(tmpdir)})

        assert status == 0

