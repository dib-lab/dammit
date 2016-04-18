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


