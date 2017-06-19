#!/usr/bin/env python
from __future__ import print_function

import json
import os
from pprint import pprint
import sys

from doit.dependency import Dependency, DbmDB

from utils import touch, datadir
from utils import run_task, run_tasks, check_status
from dammit.tasks.transdecoder import (TransDecoderLongOrfsTask,
                                       TransDecoderPredictTask)
from dammit.meta import get_config

import pytest


class TestTransDecoderTasks():

    @classmethod
    def setup_class(cls):
        cfg, _ = get_config()
        cls.longorfs_cfg = cfg['transdecoder']['longorfs']
        cls.predict_cfg = cfg['transdecoder']['predict']
        cls.extensions = ['.bed', '.pep', '.gff3', '.cds']

    def test_longorfs_task(self, tmpdir, datadir):
        with tmpdir.as_cwd():
            transcript = datadir('test-transcript.fa')
            exp_orf = datadir('test-transcript-orf.pep')

            task = TransDecoderLongOrfsTask().task(transcript,
                                                   params=self.longorfs_cfg)
            run_tasks([task], ['run'])
            output_dir = transcript + '.transdecoder_dir'

            exp_pep = open(exp_orf).read()

            pep_fn = os.path.join(output_dir, 'longest_orfs.pep')
            assert os.path.isfile(pep_fn)
            pep = open(pep_fn).read()

            assert exp_pep in pep

    def test_predict_task(self, tmpdir, datadir):
        with tmpdir.as_cwd():
            transcript = datadir('test-transcript.fa')
            pfam = datadir('test-protein-x-pfam-a.tbl')

            orf_task = TransDecoderLongOrfsTask().task(transcript,
                                                       params=self.longorfs_cfg)
            pred_task = TransDecoderPredictTask().task(transcript,
                                                       pfam,
                                                       params=self.predict_cfg)
            run_tasks([orf_task, pred_task], ['run'])
            
            pprint(tmpdir.listdir())
            for ext in self.extensions:
                fn = transcript+'.transdecoder'+ext
                assert os.path.isfile(fn)
                contents = open(fn).read()
                if ext == '.gff3':
                    assert 'mRNA' in contents
                    assert 'gene' in contents
                    assert 'CDS' in contents
                    assert 'three_prime_UTR' in contents
                    assert 'exon' in contents


