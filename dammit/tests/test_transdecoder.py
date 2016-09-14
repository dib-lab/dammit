#!/usr/bin/env python
from __future__ import print_function

import json
import os
import sys

from unittest import TestCase
from doit.dependency import Dependency, DbmDB

from utils import TemporaryDirectory, Move, TestData, touch, TemporaryFile
from utils import run_task, run_tasks, check_status
from dammit.tasks.transdecoder import (get_transdecoder_orf_task,
                                       get_transdecoder_predict_task)
from dammit.meta import get_config


class TestTransDecoderTasks(TestCase):

    @classmethod
    def setup_class(cls):
        cfg, _ = get_config()
        cls.longorfs_cfg = cfg['transdecoder']['longorfs']
        cls.predict_cfg = cfg['transdecoder']['predict']
        cls.extensions = ['.bed', '.pep', '.gff3', '.mRNA', '.cds']

    def test_longorfs_task(self):
        with TemporaryDirectory() as td:
            with Move(td):
                with TestData('test-transcript.fa', td) as transcript, \
                     TestData('test-transcript-orf.pep', td) as exp_orf:

                    task = get_transdecoder_orf_task(transcript,
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

                    orf_task = get_transdecoder_orf_task(transcript,
                                                               self.longorfs_cfg)
                    pred_task = get_transdecoder_predict_task(transcript,
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


