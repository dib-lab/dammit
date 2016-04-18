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


class TestTranscriptomeStatsTask(TestCase):

    def test_output(self):
        with TemporaryDirectory() as td:
            with Move(td):
                with TestData('test-transcript.fa', td) as transcript:

                    output_fn = os.path.join(td, 'test')
                    tsk = tasks.get_transcriptome_stats_task(transcript,
                                                             output_fn)
                    run_tasks([tsk], ['run'])

                    with open(output_fn) as fp:
                        results = json.load(fp)

                    self.assertIn('n_ambiguous', results)
                    self.assertEquals(results['n_ambiguous'], 0)


                    self.assertIn('N', results)
                    self.assertEquals(results['N'], 1)

    def test_non_acgt(self):
        with TemporaryDirectory() as td:
            with Move(td):
                with TestData('non-actg-transcripts.fa', td) as transcript:

                    output_fn = os.path.join(td, 'test')
                    tsk = tasks.get_transcriptome_stats_task(transcript,
                                                             output_fn)
                    stat = run_task(tsk)

                    self.assertEquals(stat, 2)

    def test_ambiguous_transcript(self):
        with TemporaryDirectory() as td:
            with Move(td):
                with TestData('test-transcript-N.fa', td) as transcript:

                    output_fn = os.path.join(td, 'test')
                    tsk = tasks.get_transcriptome_stats_task(transcript,
                                                             output_fn)
                    stat = run_task(tsk)

                    self.assertEquals(stat, 0)

                    print(os.listdir(td))
                    with open(output_fn) as fp:
                        results = json.load(fp)

                    self.assertIn('n_ambiguous', results)
                    self.assertEquals(results['n_ambiguous'], 1)


                    self.assertIn('N', results)
                    self.assertEquals(results['N'], 1)
                    print(results)
