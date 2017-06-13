#!/usr/bin/env python
from __future__ import print_function

import os
import stat

from dammit.app import DammitApp

from utils import datadir, runscript


PATH_BACKUP = os.environ['PATH']

def run(args, **kwargs):
    return runscript('dammit', args, **kwargs)


def test_dammit_version():
    '''Test the dammit --version command.
    '''

    from dammit.meta import __version__
    status, out, err = run(['--version'])
    print(status, out, err)
    assert status == 0
    assert out.strip() ==  'dammit {0}'.format(__version__)
