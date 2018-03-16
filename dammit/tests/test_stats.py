# Copyright (C) 2015-2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the BSD license.  See the LICENSE file for details.

import json
import os
import sys

from doit.dependency import Dependency, DbmDB

from utils import touch, datadir
from utils import run_task, run_tasks, check_status
from dammit.tasks.fastx import get_transcriptome_stats_task

class TestTranscriptomeStatsTask():

    def test_output(self, tmpdir, datadir):
        with tmpdir.as_cwd():
            transcript = datadir('test-transcript.fa')

            output_fn = str(tmpdir.join('test'))
            tsk = get_transcriptome_stats_task(transcript,
                                               output_fn)
            run_tasks([tsk], ['run'])

            with open(output_fn) as fp:
                results = json.load(fp)

            assert 'n_ambiguous' in results
            assert results['n_ambiguous'] == 0


            assert 'N' in results
            assert results['N'] == 1

    def test_non_acgt(self, tmpdir, datadir):
        with tmpdir.as_cwd():
            transcript = datadir('non-actg-transcripts.fa')
            output_fn = str(tmpdir.join('test'))

            tsk = get_transcriptome_stats_task(transcript,
                                               output_fn)
            stat = run_task(tsk)

            assert stat == 2

    def test_ambiguous_transcript(self, tmpdir, datadir):
        with tmpdir.as_cwd():
            transcript = datadir('test-transcript-N.fa')
            output_fn = str(tmpdir.join('test'))

            tsk = get_transcriptome_stats_task(transcript,
                                               output_fn)
            stat = run_task(tsk)

            with open(output_fn) as fp:
                results = json.load(fp)

            assert 'n_ambiguous' in results
            assert results['n_ambiguous'] == 1

            assert 'N' in results
            assert results['N'] == 1
            print(results)
