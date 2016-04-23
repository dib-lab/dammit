#!/usr/bin/env python
from __future__ import print_function

import json
import os
import sys

from nose.plugins.attrib import attr
from unittest import TestCase
import pandas as pd

from utils import TemporaryDirectory, Move, TestData, touch, TemporaryFile
from utils import run_task, check_status
from dammit import common
from dammit import tasks
from dammit.common import run_tasks
from dammit.infernal import get_cmpress_task as cmpress_task
from dammit.infernal import get_cmscan_task as cmscan_task
from dammit.parsers import cmscan_to_df_iter


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
                    task = cmpress_task(tf, self.cmpress_cfg)
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

                    task = cmpress_task(tf, self.cmpress_cfg)
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
                        
                    db_task = cmpress_task(cm, self.cmpress_cfg)
                    aln_task = cmscan_task(transcript, out, cm, 1.0, 1,
                                                      self.cmscan_cfg)
                    run_tasks([db_task, aln_task], ['run'])
                    print(os.listdir(td), file=sys.stderr)
                    aln = ''.join(open(out).readlines())
                    print(aln, file=sys.stderr)

                    # TODO: better correctness check
                    self.assertEquals(aln.count('accession'), 2)
                    self.assertIn('E-value', aln)

    @attr('long')
    def test_cmscan_task_multithreaded(self):
        with TemporaryDirectory() as td:
            with Move(td):
                with TestData('rnaseP-bsu.fa', td) as transcript, \
                     TestData('rnaseP-eubact.c.cm', td) as cm, \
                     TemporaryFile(td) as out_single,\
                     TemporaryFile(td) as out_multi:

                    for n_threads in (2,3,4,5):
                            
                        db_task = cmpress_task(cm, self.cmpress_cfg)
                        aln_task_single = cmscan_task(transcript, out_single, 
                                                                cm, 1.0, 1,
                                                                self.cmscan_cfg)
                        aln_task_multi = cmscan_task(transcript, out_multi, 
                                                                cm, 1.0,
                                                                n_threads,
                                                                self.cmscan_cfg)
                        run_tasks([db_task, aln_task_single], ['run'])
                        run_task(aln_task_multi)

                        alns_single = pd.concat(cmscan_to_df_iter(out_single))
                        alns_multi = pd.concat(cmscan_to_df_iter(out_multi))

                        self.assertTrue(all(alns_single['e_value'].sort_values() == \
                                            alns_multi['e_value'].sort_values()))


