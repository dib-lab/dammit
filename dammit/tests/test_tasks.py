#!/usr/bin/env python
from __future__ import print_function

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

Shit your 'scripts ain't differential --
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
                    print('PATH:', os.environ['PATH'], file=sys.stderr)

                    
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
                    aln_task = tasks.get_lastal_task(tr, prot, out, True, 1,
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
                    aln_task = tasks.get_lastal_task(prot, prot, out, False, 1,
                                                     self.lastal_cfg)
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
                    aln_task = tasks.get_lastal_task(prot, prot, out, False, 1,
                                                     self.lastal_cfg)
                    # Run it once
                    run_tasks([db_task, aln_task], ['run'])

                    # Now run again and check the status
                    run_tasks([aln_task], ['run'])
                    status = check_status(aln_task)
                    self.assertEquals(status.status, 'up-to-date')
