import os
import pytest
from utils import run, runscript


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
        status, out, err = snakemake_rule('busco/busco.rule',
                                           target='run_busco',
                                           config={'DATA_DIR': str(tmpdir)},
                                           extra_args=["-n"])

        assert status == 0


#@pytest.mark.long
#def test_busco(snakemake_rule, tmpdir, datadir):
#    with tmpdir.as_cwd():
#        datadir('target.fa')
#        status, out, err = snakemake_rule('busco/busco.rule',
#                                           target='run_busco',
#                                           config={'DATA_DIR': str(tmpdir)})
#
#        print(out)
#        assert os.path.isfile('txome_busco/full_table_txome_busco.tsv')
#        assert os.path.isfile('test-busco.log')
#        print(open('test-busco.log').read())
#        #assert '[ok]' in open('test-busco.log').read()

