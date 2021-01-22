import os
import shutil
import sys

from dammit.utils import create_dirs
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
        item._nodeid = '{0}: {1} {2}'.format(item.name, pref, suf)


@pytest.fixture(scope="session")
def monkeysession(request):
    from _pytest.monkeypatch import MonkeyPatch
    mpatch = MonkeyPatch()
    yield mpatch
    mpatch.undo()


@pytest.fixture(scope='session', autouse=True)
def setup_test_environment(monkeysession, tmpdir_factory):
    testing_temp_base_dir = os.getenv('DAMMIT_TESTING_TEMP_BASE_DIR')
    if testing_temp_base_dir is None:
        testing_temp_base_dir = tmpdir_factory.mktemp('dammit-testing-temp-base')
    print('Testing temp base directory:', testing_temp_base_dir, file=sys.stderr)
    testing_temp_temp_dir = os.path.join(testing_temp_base_dir, 'temp')
    testing_conda_dir     = os.path.join(testing_temp_base_dir, 'envs')

    create_dirs([testing_temp_temp_dir, testing_conda_dir])
    monkeysession.setenv('DAMMIT_TEMP_DIR', str(testing_temp_temp_dir))
    monkeysession.setenv('DAMMIT_CONDA_DIR', str(testing_conda_dir))
    return testing_temp_base_dir


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
def snakemake_rule(setup_test_environment):
    conda_env_dir = os.path.join(setup_test_environment, 'envs')

    def run(rule_path, target=None, config=None, n_threads=1, extra_args = [], **kwargs):
        if '--config' in extra_args:
            raise RuntimeError('pass --config to config keyword')
        rule_path = os.path.join(__wrappers__, rule_path)

        args = []
        if config is not None:
            args.append('--config')
            args.append(' '.join((f'{k}={v}' for k, v in config.items())))

        args.extend(['-s', rule_path, '--use-conda', '--conda-prefix', conda_env_dir, '-p', '-j', str(n_threads)])
        args.extend(extra_args)
        if target is not None:
            if isinstance(target, str):
                args.append(target)
            else:
                args.extend(target)
        print(' '.join(args))
        return runscript('snakemake', args, **kwargs)
    return run
