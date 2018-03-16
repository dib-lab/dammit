# Copyright (C) 2015-2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the BSD license.  See the LICENSE file for details.

import json
import os
import sys

import pandas as pd

from utils import touch, datadir
from utils import run_task, run_tasks, check_status
from dammit.tasks.infernal import CMPressTask, CMScanTask
from dammit.fileio.infernal import InfernalParser
from dammit.meta import get_config

import pytest

class TestInfernalTasks():

    @classmethod
    def setup_class(cls):
        cfg, _ = get_config()
        cls.cmscan_cfg = cfg['infernal']['cmscan']
        cls.cmpress_cfg = cfg['infernal']['cmscan']
        cls.extensions = ['.i1f', '.i1i', '.i1m', '.i1p']

    def test_cmpress_task(self, tmpdir, datadir):
        with tmpdir.as_cwd():
            tf = datadir('test-covariance-model.cm')
            task = CMPressTask().task(tf, params=self.cmpress_cfg)
            run_tasks([task], ['run'])
            status = check_status(task)
            print(os.listdir(), file=sys.stderr)
            
            for ext in self.extensions:
                assert os.path.isfile(tf + ext)

            assert status.status == 'up-to-date'

    def test_cmpress_task_existing(self, tmpdir, datadir):
        with tmpdir.as_cwd():
            tf = datadir('test-covariance-model.cm')
            for ext in self.extensions:
                touch(tf + ext)

            task = CMPressTask().task(tf, params=self.cmpress_cfg)
            run_tasks([task], ['run'])
            print(os.listdir(), file=sys.stderr)
            print(task, file=sys.stderr)
            status = check_status(task)
            
            assert status.status == 'up-to-date'

    def test_cmscan_task(self, tmpdir, datadir):
        with tmpdir.as_cwd():
            transcript = datadir('test-transcript.fa')
            cm = datadir('test-covariance-model.cm')
            out = str(tmpdir.join('test.out'))
                
            db_task = CMPressTask().task(cm, params=self.cmpress_cfg)
            aln_task = CMScanTask().task(transcript, out, cm, 
                                         cutoff=1.0, n_threads=1)
            run_tasks([db_task, aln_task], ['run'])
            print(os.listdir(), file=sys.stderr)
            aln = ''.join(open(out).readlines())
            print(aln, file=sys.stderr)

            # TODO: better correctness check
            assert aln.count('accession') == 2
            assert 'E-value' in aln

    def test_cmscan_task_multithreaded(self, tmpdir, datadir):
        with tmpdir.as_cwd():
            transcript = datadir('rnaseP-bsu.fa')
            cm = datadir('rnaseP-eubact.c.cm')
            out_single = str(tmpdir.join('single'))
            out_multi = str(tmpdir.join('multi'))

            for n_threads in (2,3,4,5):
                    
                db_task = CMPressTask().task(cm, params=self.cmpress_cfg)
                aln_tasks_single = CMScanTask().task(transcript, out_single, 
                                                     cm, cutoff=1.0, 
                                                     n_threads=1)
                aln_tasks_multi = CMScanTask().task(transcript, out_multi, 
                                                    cm, cutoff=1.0,
                                                    n_threads=n_threads)
                run_tasks([db_task, aln_tasks_single], ['run'])
                run_task(aln_tasks_multi)

                alns_single = pd.concat(InfernalParser(out_single))
                alns_multi = pd.concat(InfernalParser(out_multi))

                assert all(alns_single['e_value'].sort_values() == \
                           alns_multi['e_value'].sort_values())


