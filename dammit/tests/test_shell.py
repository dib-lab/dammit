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
