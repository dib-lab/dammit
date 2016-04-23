#!/usr/bin/env python
from __future__ import print_function

import json
import os
import sys

from unittest import TestCase


from utils import TemporaryDirectory, Move, TestData, touch, TemporaryFile
from utils import run_task, check_status
from dammit import common
from dammit import tasks
from dammit.last import get_lastdb_task as lastdb_task
from dammit.last import get_lastal_task as lastal_task
from dammit.common import run_tasks
from dammit.fileio.maf import MafParser


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
                    task = lastdb_task(tf, tf, self.lastdb_cfg,
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
                    task = lastdb_task(tf, tf, self.lastdb_cfg,
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

                    task = lastdb_task(tf, tf, self.lastdb_cfg,
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
                    db_task = lastdb_task(prot, prot, self.lastdb_cfg)
                    aln_task = lastal_task(tr, prot, out, self.lastal_cfg, 
                                                     translate=True, 
                                                     cutoff=None)
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
                    db_task = lastdb_task(prot, prot, self.lastdb_cfg)
                    aln_task = lastal_task(prot, prot, out,
                                                     self.lastal_cfg,
                                                     translate=False,
                                                     cutoff=None)
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

    def test_lastal_task_multithreaded(self):
        with TemporaryDirectory() as td:
            with Move(td):
                with TestData('test-protein.fa', td) as prot, \
                     TestData('pom.50.fa', td) as tr, \
                     TemporaryFile(td) as out_single,\
                     TemporaryFile(td) as out_multi:
                        
                    for n_threads in (3,4,5):
                        print(os.listdir(td), file=sys.stderr)

                        db_task = lastdb_task(prot, prot, self.lastdb_cfg)
                        aln_task_single = lastal_task(tr, prot, out_single, 
                                                                self.lastal_cfg, 
                                                                translate=True, 
                                                                cutoff=None)

                        aln_task_multi = lastal_task(tr, prot, out_multi, 
                                                                self.lastal_cfg, 
                                                                translate=True, 
                                                                cutoff=None,
                                                                n_threads=n_threads)
                        run_tasks([db_task, aln_task_single, aln_task_multi], 
                                  ['run'])

                        alns_single = MafParser(out_single).read()
                        alns_multi = MafParser(out_multi).read()

                        self.assertTrue(all(alns_single['E'].sort_values() == \
                                        alns_multi['E'].sort_values()))

    def test_lastal_task_uptodate(self):
        with TemporaryDirectory() as td:
            with Move(td):
                with TestData('test-protein.fa', td) as prot, \
                     TemporaryFile(td) as out:
                        
                    print(os.listdir(td), file=sys.stderr)
                    db_task = lastdb_task(prot, prot, self.lastdb_cfg)
                    aln_task = lastal_task(prot, prot, out,
                                                     self.lastal_cfg,
                                                     translate=False,
                                                     cutoff=None)
                    # Run it once
                    run_tasks([db_task, aln_task], ['run'])

                    # Now run again and check the status
                    run_tasks([aln_task], ['run'])
                    status = check_status(aln_task)
                    self.assertEquals(status.status, 'up-to-date')
