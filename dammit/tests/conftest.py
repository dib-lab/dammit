import os
import shutil

from dammit.config import DEFAULT_TEMP_DIR
from dammit.meta import __wrappers__
from .utils import runscript

import pytest


def pytest_itemcollected(item):
    par = item.parent.obj
    node = item.obj
    pref = par.__doc__.strip() if par.__doc__ else par.__class__.__name__
    if pref == 'module':
        pref = ''
    suf = node.__doc__.strip() if node.__doc__ else node.__name__
    if pref or suf:
        item._nodeid = ' '.join((pref, suf))


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
    def run(rule_path, target=None, config=None, extra_args = [], **kwargs):
        if '--config' in extra_args:
            raise RuntimeError('pass --config to config keyword')
        rule_path = os.path.join(__wrappers__, rule_path)

        args = []
        if config is not None:
            args.append('--config')
            args.append(' '.join((f'{k}={v}' for k, v in config.items())))

        args.extend(['-s', rule_path, '--use-conda', '--conda-prefix', conda_env_dir, '-j', '1'])
        args.extend(extra_args)
        if target is not None:
            if isinstance(target, str):
                args.append(target)
            else:
                args.extend(target)
        print(' '.join(args))
        return runscript('snakemake', args, **kwargs)
    return run
