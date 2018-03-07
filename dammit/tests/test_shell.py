# Copyright (C) 2015-2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the BSD license.  See the LICENSE file for details.

import os
import stat
import pandas as pd
import pytest

from dammit.tasks.shell import get_download_task, hashfile

from utils import datadir, runscript

def test_hashfile_filenotfounderror(tmpdir):
    bad_file_path = str(tmpdir.join('not_a_file'))
    with pytest.raises(RuntimeError):
        hashfile(bad_file_path)
