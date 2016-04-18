#!/usr/bin/env python
from __future__ import print_function

import json
import os
import sys

from unittest import TestCase
from doit.dependency import Dependency, DbmDB

from utils import TemporaryDirectory, Move, TestData, touch, TemporaryFile
from utils import run_task, check_status
from dammit import common
from dammit import tasks
from dammit.common import run_tasks

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


