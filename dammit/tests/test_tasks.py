#!/usr/bin/env python
from __future__ import print_function

import os
import sys

from unittest import TestCase
from doit.dependency import Dependency, DbmDB

from utils import TemporaryDirectory, Move, TestData
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
        cls.lastdb_cfg = common.CONFIG['settings']['last']['lastal']

    def test_lastdb_nucl_task(self):
        status = None
        with TemporaryDirectory() as td:
            with Move(td):
                with TestData('test-transcript.fa', td) as tf:
                    task = tasks.get_lastdb_task(tf, tf, self.lastdb_cfg,
                                                 prot=False)
                    run_tasks([task], ['run'])
                    status = check_status(task)
                    print(os.listdir(td), file=sys.stderr)
                    print('PATH:', os.environ['PATH'], file=sys.stderr)

                    
                    for ext in ['.bck', '.des', '.prj', '.sds', 
                                '.ssp', '.suf', '.tis']:
                        self.assertTrue(os.path.isfile(tf + ext))

                    self.assertEquals(status.status, 'up-to-date')
