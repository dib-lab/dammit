# Copyright (C) 2015-2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the BSD license.  See the LICENSE file for details.

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
