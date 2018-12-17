# Copyright (C) 2015-2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the BSD license.  See the LICENSE file for details.

import json
import os
import sys

import pandas as pd

from utils import datadir, touch
from utils import run_task, run_tasks, check_status

from dammit.tasks.hmmer import HMMScanTask, HMMPressTask
from dammit.fileio.hmmer import HMMerParser
from dammit.meta import get_config

import pytest

class TestHMMERTasks():

    @classmethod
    def setup_class(cls):
        cfg, _ = get_config()
        cls.hmmscan_cfg = cfg['hmmer']['hmmscan']
        cls.hmmpress_cfg = cfg['hmmer']['hmmpress']
        cls.extensions = ['.h3f', '.h3i', '.h3m', '.h3p']

    def test_hmmpress_task(self, tmpdir, datadir):
        with tmpdir.as_cwd():
            tf = datadir('test-profile.hmm')
            task = HMMPressTask().task(tf)
            run_tasks([task], ['run'])
            status = check_status(task)
            print(os.listdir(), file=sys.stderr)
            
            for ext in self.extensions:
                assert os.path.isfile(tf + ext)

            assert status.status == 'up-to-date'

    def test_hmmpress_task_existing(self, tmpdir, datadir):
        with tmpdir.as_cwd():        
            tf = datadir('test-profile.hmm')
            for ext in self.extensions:
                touch(tf + ext)
            task = HMMPressTask().task(tf)
            run_tasks([task], ['run'])
            print(os.listdir(), file=sys.stderr)
            print(task, file=sys.stderr)
            status = check_status(task)
            
            assert status.status == 'up-to-date'

    def test_hmmscan_task(self, tmpdir, datadir):
        with tmpdir.as_cwd():
            prot = datadir('test-protein.fa')
            hmm = datadir('test-profile.hmm')
            out = str(tmpdir.join('test.out'))
                        
            db_task = HMMPressTask().task(hmm, params=self.hmmpress_cfg)
            aln_task = HMMScanTask().task(prot, out, hmm, 
                                          cutoff=1.0, n_threads=1)

            run_tasks([db_task, aln_task], ['run'])
            print(os.listdir(), file=sys.stderr)
            aln = open(out).read()
            print(aln)

            assert aln.count('accession') == 2
            assert 'i-Evalue' in aln

    def test_hmmscan_task_multithreaded(self, tmpdir, datadir):
        with tmpdir.as_cwd():
            prot = datadir('20aa-alitest.fa')
            hmm = datadir('20aa.hmm')
            out_single = str(tmpdir.join('out-single'))
            out_multi = str(tmpdir.join('out-multi'))
    
            for n_threads in (2,3,4,5):
                db_task = HMMPressTask().task(hmm, params=self.hmmpress_cfg)
                aln_task_single = HMMScanTask().task(prot, out_single, 
                                                     hmm, cutoff=1.0, 
                                                     n_threads=1)
                aln_task_multi = HMMScanTask().task(prot, out_multi,
                                                    hmm, cutoff=1.0, 
                                                    n_threads=n_threads)
                run_tasks([db_task, aln_task_single], ['run'])
                run_tasks([aln_task_multi], ['run'])
                print(os.listdir(), file=sys.stderr)

                print(open(out_single).read())
                alns_single = pd.concat(HMMerParser(out_single))
                alns_multi = pd.concat(HMMerParser(out_multi))

                assert all(alns_single['domain_i_evalue'].sort_values().reset_index(drop=True) == \
                           alns_multi['domain_i_evalue'].sort_values().reset_index(drop=True))


