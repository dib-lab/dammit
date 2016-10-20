from unittest import TestCase

from nose.plugins.attrib import attr
import os
import stat
import pandas as pd

from dammit.app import DammitApp
from dammit import dependencies

from utils import TemporaryDirectory, TestData, runscript

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

names = ['TransDecoder',
         'crb-blast',
         'LAST',
         'BUSCO',
         'BLAST+',
         'Infernal',
         'HMMER']

execs = ['hmmscan',
         'hmmpress',
         'cmscan',
         'cmpress',
         'blastp',
         'blastx',
         'tblastn',
         'makeblastdb',
         'BUSCO_v1.1b1.py',
         'TransDecoder.LongOrfs',
         'TransDecoder.Predict',
         'lastdb',
         'crb-blast']

PATH_BACKUP = os.environ['PATH']

def run(args, **kwargs):
    return runscript('dammit', args, **kwargs)


def mock_last(directory):
    fn = os.path.join(directory, 'lastal')
    with open(fn, 'w') as fp:
        fp.write('#!/bin/sh\n\necho "lastal 611"')
    os.chmod(fn, os.stat(fn).st_mode | stat.S_IEXEC)
    print('added mock lastal to', fn)


class TestDependencies(TestCase):

    @classmethod
    def setup_class(cls):
        cls.maxDiff = None
        #cls.PATH_BACKUP = os.environ['PATH']
        #print('Setup. PATH backup:', cls.PATH_BACKUP)

    @classmethod
    def teardown_class(cls):
        os.environ['PATH'] = PATH_BACKUP

    @staticmethod
    def add_execs_to_path(tempdir):
        os.environ['PATH'] = ''
        for exe in execs:
            fn = os.path.join(tempdir, exe)
            open(fn, 'a').close()
            os.chmod(fn, stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
        os.environ['PATH'] += os.pathsep + tempdir
        mock_last(tempdir)
        print('Added execs:', os.environ['PATH'])
        print('dir:', os.listdir(tempdir))

    def test_get_all_statuses_default_nodeps(self):

        os.environ['PATH'] = ''
        handler = dependencies.get_handler()
        is_fulfilled, unfulfilled = handler.get_all_statuses()
        self.assertFalse(is_fulfilled)
        for name, msg in unfulfilled.items():
            self.assertIn(name, names)

    def test_get_all_statuses_default_alldeps(self):
        handler = dependencies.get_handler()
        with TemporaryDirectory() as tempdir:
            TestDependencies.add_execs_to_path(tempdir)

            is_fulfilled, unfulfilled = handler.get_all_statuses()
            print(unfulfilled)
            self.assertTrue(is_fulfilled)
            for name, msg in unfulfilled.items():
                self.assertIn(name, names)

    def test_print_all_statuses_default_no_deps(self):
        os.environ['PATH'] = ''
        handler = dependencies.get_handler()
        out = StringIO()

        is_fulfilled, unfulfilled = handler.print_all_statuses(out=out)
        output = out.getvalue()

        self.assertFalse(is_fulfilled)
        self.assertIn('Some dependencies unfulfilled', output)

    def test_print_all_statuses_default_alldeps(self):
        handler = dependencies.get_handler()
        with TemporaryDirectory() as tempdir:
            TestDependencies.add_execs_to_path(tempdir)
            out = StringIO()

            is_fulfilled, unfulfilled = handler.print_all_statuses(out=out)
            print(unfulfilled)
            output = out.getvalue()

            self.assertTrue(is_fulfilled)
            self.assertIn('All dependencies fulfilled!', output)

    def test_check_or_fail_default_no_deps(self):
        os.environ['PATH'] = ''
        handler = dependencies.get_handler()
        out = StringIO()
        with self.assertRaises(SystemExit):
            handler.check_or_fail(out=out)
        output = out.getvalue()
        self.assertIn('Must install', output)

    def test_check_or_fail_default_alldeps(self):
        handler = dependencies.get_handler()
        with TemporaryDirectory() as tempdir:
            TestDependencies.add_execs_to_path(tempdir)
            out = StringIO()
            
            try:
                handler.check_or_fail(out=out)
            except SystemExit:
                assert False, 'check should have succeeded'

class TestDammitDependencies(TestCase):

    @classmethod
    def teardown_class(cls):
        os.environ['PATH'] = PATH_BACKUP

    def test_dammit_dependencies(self):
        '''Test the dependencies subcommand.
        '''
        print(os.environ['PATH'])
        status, out, err = run(['dependencies'])
        self.assertEquals(status, 0)

