import os

from dammit.config import DEFAULT_TEMP_DIR
from dammit.meta import __wrappers__
from utils import runscript

import pytest


@pytest.fixture
def conda_env_dir():
    TEST_ENV_DIR = os.path.join(DEFAULT_TEMP_DIR, 'test-envs')
    os.makedirs(TEST_ENV_DIR, exist_ok=True)
    return TEST_ENV_DIR


@pytest.fixture
def datadir(tmpdir, request):
    '''
    Fixture responsible for locating the test data directory and copying it
    into a temporary directory.
    '''
    filename   = request.module.__file__
    test_dir   = os.path.dirname(filename)
    data_dir   = os.path.join(test_dir, 'test-data')

    def getter(filename, as_str=True):
        filepath = tmpdir.join(filename)
        shutil.copyfile(os.path.join(data_dir, filename),
                        filepath)
        if as_str:
            return str(filepath)
        return filepath

    return getter


@pytest.fixture
def snakemake_rule(conda_env_dir):
    def run(rule_path, target=None, extra_args = [], **kwargs):
        rule_path = os.path.join(__wrappers__, rule_path)
        args = ['-s', rule_path, '--use-conda', '--conda-prefix', conda_env_dir, '-j', '1']
        if target is not None:
            if isinstance(target, str):
                args.append(target)
            else:
                args.extend(target)
        args.extend(extra_args)
        return runscript('snakemake', args, **kwargs)
    return run