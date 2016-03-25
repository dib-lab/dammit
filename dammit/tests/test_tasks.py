#!/usr/bin/env python
from __future__ import print_function

import json
import os
import sys

from unittest import TestCase
from doit.dependency import Dependency, DbmDB

from utils import TemporaryDirectory, Move, TestData, touch, TemporaryFile
from dammit import common
from dammit import tasks
from dammit.common import run_tasks

'''
BATCH EFFECTS -- The Notorious A.T.G.
-------------------------------------

Shit your 'scripts ain't differential -
they're preferential
(-ly selected!)
from you read seq
to your DEseq
your biases are your analyses fallacies
your matrices
are make-believe
your shallow e-values swallowing your logic and your counts bouncin' --
you got BATCH EFFECTS
(*batch effects*)
BATCH EFFECTS
(*batch effects*)

'''


def check_status(task, dep_file='.doit.db'):
    mgr = Dependency(DbmDB, os.path.abspath(dep_file))
    status = mgr.get_status(task, [task])
    return status


def run_task(task, cmd='run', verbosity=2):
    return run_tasks([task], [cmd], config={'verbosity': verbosity})


class TestLASTTasks(TestCase):

    @classmethod
    def setup_class(cls):
        cls.lastdb_cfg = common.CONFIG['settings']['last']['lastdb']
        cls.lastal_cfg = common.CONFIG['settings']['last']['lastal']
        cls.extensions = ['.bck', '.des', '.prj', '.sds', '.ssp', '.suf', '.tis']

    def test_lastdb_task_nucl(self):
        with TemporaryDirectory() as td:
            with Move(td):
                with TestData('test-transcript.fa', td) as tf:
                    task = tasks.get_lastdb_task(tf, tf, self.lastdb_cfg,
                                                 prot=False)
                    run_tasks([task], ['run'])
                    status = check_status(task)
                    print(os.listdir(td), file=sys.stderr)
                    print('PATH:', os.environ['PATH'], file=sys.stderr)

                    
                    for ext in self.extensions:
                        self.assertTrue(os.path.isfile(tf + ext))

                    self.assertEquals(status.status, 'up-to-date')

    def test_lastdb_task_prot(self):
        with TemporaryDirectory() as td:
            with Move(td):
                with TestData('test-protein.fa', td) as tf:
                    task = tasks.get_lastdb_task(tf, tf, self.lastdb_cfg,
                                                 prot=True)
                    run_tasks([task], ['run'])
                    status = check_status(task)
                    print(os.listdir(td), file=sys.stderr)
                    
                    for ext in self.extensions:
                        self.assertTrue(os.path.isfile(tf + ext))

                    self.assertEquals(status.status, 'up-to-date')

    def test_lastdb_task_existing(self):
        with TemporaryDirectory() as td:
            with Move(td):
                with TestData('test-protein.fa', td) as tf:
                    for ext in self.extensions:
                        touch(tf + ext)

                    task = tasks.get_lastdb_task(tf, tf, self.lastdb_cfg,
                                                 prot=True)
                    run_tasks([task], ['run'])
                    print(os.listdir(td), file=sys.stderr)
                    print(task, file=sys.stderr)
                    status = check_status(task)

                    self.assertEquals(status.status, 'up-to-date')

    def test_lastal_task_nucl_x_prot(self):
        with TemporaryDirectory() as td:
            with Move(td):
                with TestData('test-protein.fa', td) as prot, \
                     TestData('test-transcript.fa', td) as tr, \
                     TemporaryFile(td) as out:
                        
                    print(os.listdir(td), file=sys.stderr)
                    db_task = tasks.get_lastdb_task(prot, prot, self.lastdb_cfg)
                    aln_task = tasks.get_lastal_task(tr, prot, out, True, None, 1,
                                                     self.lastal_cfg)
                    run_tasks([db_task, aln_task], ['run'])

                    aln = ''.join(open(out).readlines())
                    print(aln, file=sys.stderr)

                    self.assertIn('SPAC212_RecQ_type_DNA_helicase_PROTEIN', 
                                  aln)
                    self.assertIn('SPAC212_RecQ_type_DNA_helicase_TRANSCRIPT',
                                  aln)
                    self.assertIn('lambda', aln, 
                                  msg='lambda missing, wrong LAST version?')
                    

    def test_lastal_task_prot_x_prot(self):
        with TemporaryDirectory() as td:
            with Move(td):
                with TestData('test-protein.fa', td) as prot, \
                     TemporaryFile(td) as out:
                        
                    print(os.listdir(td), file=sys.stderr)
                    db_task = tasks.get_lastdb_task(prot, prot, self.lastdb_cfg)
                    aln_task = tasks.get_lastal_task(prot, prot, out, False,
                                                     None, 1, self.lastal_cfg)
                    run_tasks([db_task, aln_task], ['run'])

                    aln = ''.join(open(out).readlines())
                    print(aln, file=sys.stderr)

                    self.assertEquals(
                            aln.count('SPAC212_RecQ_type_DNA_helicase_PROTEIN'),
                            2)
                    self.assertIn('EG2=0', aln)
                    self.assertIn('E=0', aln)
                    self.assertIn('lambda', aln, 
                                  msg='lambda missing, wrong LAST version?')

    def test_lastal_task_uptodate(self):
        with TemporaryDirectory() as td:
            with Move(td):
                with TestData('test-protein.fa', td) as prot, \
                     TemporaryFile(td) as out:
                        
                    print(os.listdir(td), file=sys.stderr)
                    db_task = tasks.get_lastdb_task(prot, prot, self.lastdb_cfg)
                    aln_task = tasks.get_lastal_task(prot, prot, out, False,
                                                     None, 1, self.lastal_cfg)
                    # Run it once
                    run_tasks([db_task, aln_task], ['run'])

                    # Now run again and check the status
                    run_tasks([aln_task], ['run'])
                    status = check_status(aln_task)
                    self.assertEquals(status.status, 'up-to-date')


class TestHMMERTasks(TestCase):

    @classmethod
    def setup_class(cls):
        cls.hmmscan_cfg = common.CONFIG['settings']['hmmer']['hmmscan']
        cls.hmmpress_cfg = common.CONFIG['settings']['hmmer']['hmmscan']
        cls.extensions = ['.h3f', '.h3i', '.h3m', '.h3p']

    def test_hmmpress_task(self):
        with TemporaryDirectory() as td:
            with Move(td):
                with TestData('test-profile.hmm', td) as tf:
                    task = tasks.get_hmmpress_task(tf, self.hmmpress_cfg)
                    run_tasks([task], ['run'])
                    status = check_status(task)
                    print(os.listdir(td), file=sys.stderr)
                    
                    for ext in self.extensions:
                        self.assertTrue(os.path.isfile(tf + ext))

                    self.assertEquals(status.status, 'up-to-date')
    
    def test_hmmpress_task_existing(self):
        with TemporaryDirectory() as td:
            with Move(td):
                with TestData('test-profile.hmm', td) as tf:
                    for ext in self.extensions:
                        touch(tf + ext)

                    task = tasks.get_hmmpress_task(tf, self.hmmpress_cfg)
                    run_tasks([task], ['run'])
                    print(os.listdir(td), file=sys.stderr)
                    print(task, file=sys.stderr)
                    status = check_status(task)
                    
                    self.assertEquals(status.status, 'up-to-date')

    def test_hmmscan_task(self):
        with TemporaryDirectory() as td:
            with Move(td):
                with TestData('test-protein.fa', td) as prot, \
                     TestData('test-profile.hmm', td) as hmm, \
                     TemporaryFile(td) as out:
                        
                    db_task = tasks.get_hmmpress_task(hmm, self.hmmpress_cfg)
                    aln_task = tasks.get_hmmscan_task(prot, out, hmm, 1.0, 1,
                                                      self.hmmscan_cfg)
                    run_tasks([db_task, aln_task], ['run'])
                    print(os.listdir(td), file=sys.stderr)
                    aln = ''.join(open(out).readlines())
                    print(aln, file=sys.stderr)

                    self.assertEquals(aln.count('accession'), 2)
                    self.assertIn('i-Evalue', aln)


class TestInfernalTasks(TestCase):

    @classmethod
    def setup_class(cls):
        cls.cmscan_cfg = common.CONFIG['settings']['infernal']['cmscan']
        cls.cmpress_cfg = common.CONFIG['settings']['infernal']['cmscan']
        cls.extensions = ['.i1f', '.i1i', '.i1m', '.i1p']

    def test_cmpress_task(self):
        with TemporaryDirectory() as td:
            with Move(td):
                with TestData('test-covariance-model.cm', td) as tf:
                    task = tasks.get_cmpress_task(tf, self.cmpress_cfg)
                    run_tasks([task], ['run'])
                    status = check_status(task)
                    print(os.listdir(td), file=sys.stderr)
                    
                    for ext in self.extensions:
                        self.assertTrue(os.path.isfile(tf + ext))

                    self.assertEquals(status.status, 'up-to-date')
    
    def test_cmpress_task_existing(self):
        with TemporaryDirectory() as td:
            with Move(td):
                with TestData('test-covariance-model.cm', td) as tf:
                    for ext in self.extensions:
                        touch(tf + ext)

                    task = tasks.get_cmpress_task(tf, self.cmpress_cfg)
                    run_tasks([task], ['run'])
                    print(os.listdir(td), file=sys.stderr)
                    print(task, file=sys.stderr)
                    status = check_status(task)
                    
                    self.assertEquals(status.status, 'up-to-date')

    def test_cmscan_task(self):
        with TemporaryDirectory() as td:
            with Move(td):
                with TestData('test-transcript.fa', td) as transcript, \
                     TestData('test-covariance-model.cm', td) as cm, \
                     TemporaryFile(td) as out:
                        
                    db_task = tasks.get_cmpress_task(cm, self.cmpress_cfg)
                    aln_task = tasks.get_cmscan_task(transcript, out, cm, 1.0, 1,
                                                      self.cmscan_cfg)
                    run_tasks([db_task, aln_task], ['run'])
                    print(os.listdir(td), file=sys.stderr)
                    aln = ''.join(open(out).readlines())
                    print(aln, file=sys.stderr)

                    # TODO: better correctness check
                    self.assertEquals(aln.count('accession'), 2)
                    self.assertIn('E-value', aln)


class TestTransDecoderTasks(TestCase):

    @classmethod
    def setup_class(cls):
        cls.longorfs_cfg = common.CONFIG['settings']['transdecoder']['longorfs']
        cls.predict_cfg = common.CONFIG['settings']['transdecoder']['predict']
        cls.extensions = ['.bed', '.pep', '.gff3', '.mRNA', '.cds']

    def test_longorfs_task(self):
        with TemporaryDirectory() as td:
            with Move(td):
                with TestData('test-transcript.fa', td) as transcript, \
                     TestData('test-transcript-orf.pep', td) as exp_orf:

                    task = tasks.get_transdecoder_orf_task(transcript,
                                                           self.longorfs_cfg)
                    run_tasks([task], ['run'])
                    output_dir = transcript + '.transdecoder_dir'

                    exp_pep = open(exp_orf).read()

                    pep_fn = os.path.join(output_dir, 'longest_orfs.pep')
                    self.assertTrue(os.path.isfile(pep_fn))
                    pep = open(pep_fn).read()

                    self.assertIn(exp_pep, pep)

    def test_predict_task(self):
        with TemporaryDirectory() as td:
            with Move(td):
                with TestData('test-transcript.fa', td) as transcript, \
                     TestData('test-protein-x-pfam-a.tbl', td) as pfam:

                    orf_task = tasks.get_transdecoder_orf_task(transcript,
                                                               self.longorfs_cfg)
                    pred_task = tasks.get_transdecoder_predict_task(transcript,
                                                                    pfam,
                                                                    self.predict_cfg)
                    run_tasks([orf_task, pred_task], ['run'])
                    
                    for ext in self.extensions:
                        fn = os.path.join(td, transcript+'.transdecoder'+ext)
                        self.assertTrue(os.path.isfile(fn))
                        contents = open(fn).read()
                        if ext == '.gff3':
                            self.assertIn('mRNA', contents)
                            self.assertIn('gene', contents)
                            self.assertIn('CDS', contents)
                            self.assertIn('three_prime_UTR', contents)
                            self.assertIn('exon', contents)


class TestTranscriptomeStatsTask(TestCase):

    def test_output(self):
        with TemporaryDirectory() as td:
            with Move(td):
                with TestData('test-transcript.fa', td) as transcript:

                    output_fn = os.path.join(td, 'test')
                    tsk = tasks.get_transcriptome_stats_task(transcript,
                                                             output_fn)
                    run_tasks([tsk], ['run'])

                    with open(output_fn) as fp:
                        results = json.load(fp)

                    self.assertIn('n_ambiguous', results)
                    self.assertEquals(results['n_ambiguous'], 0)


                    self.assertIn('N', results)
                    self.assertEquals(results['N'], 1)

    def test_non_acgt(self):
        with TemporaryDirectory() as td:
            with Move(td):
                with TestData('non-actg-transcripts.fa', td) as transcript:

                    output_fn = os.path.join(td, 'test')
                    tsk = tasks.get_transcriptome_stats_task(transcript,
                                                             output_fn)
                    stat = run_task(tsk)

                    self.assertEquals(stat, 2)

    def test_ambiguous_transcript(self):
        with TemporaryDirectory() as td:
            with Move(td):
                with TestData('test-transcript-N.fa', td) as transcript:

                    output_fn = os.path.join(td, 'test')
                    tsk = tasks.get_transcriptome_stats_task(transcript,
                                                             output_fn)
                    stat = run_task(tsk)

                    self.assertEquals(stat, 0)

                    print(os.listdir(td))
                    with open(output_fn) as fp:
                        results = json.load(fp)

                    self.assertIn('n_ambiguous', results)
                    self.assertEquals(results['n_ambiguous'], 1)


                    self.assertIn('N', results)
                    self.assertEquals(results['N'], 1)
                    print(results)
