#!/usr/bin/env python
from __future__ import print_function
from __future__ import absolute_import

import nose
import os
from io import StringIO
import logging
import traceback
import shutil
import stat
import sys
import warnings as _warnings
from pkg_resources import Requirement, resource_filename, ResolutionError
from tempfile import mkdtemp

from doit.dependency import Dependency, DbmDB

from dammit import log
log.start_logging(test=True)
logger = logging.getLogger('Tests')

'''
BATCH EFFECTS -- The Notorious A.T.G.
-------------------------------------

Shit your 'scripts ain't differential -
they're preferential
(-ly selected!)
from you read seq
to your DEseq
your biases are your analyses fallacies
your matrices
are make-believe
your shallow e-values swallowing your logic and your counts bouncin' --
you got BATCH EFFECTS
(*batch effects*)
BATCH EFFECTS
(*batch effects*)

'''


try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO


def check_status(task, tasks=None, dep_file='.doit.db'):
    if tasks is None:
        tasks = [task]
    mgr = Dependency(DbmDB, os.path.abspath(dep_file))
    status = mgr.get_status(task, tasks)
    return status


def run_tasks(tasks, args, config={'verbosity': 0}):
    
    if type(tasks) is not list:
        raise TypeError('tasks must be a list')
   
    class Loader(TaskLoader):
        @staticmethod
        def load_tasks(cmd, opt_values, pos_args):
            return tasks, config
   
    return DoitMain(Loader()).run(args)


def run_task(task, cmd='run', verbosity=2):
    return run_tasks([task], [cmd], config={'verbosity': verbosity})


def touch(filename):
    '''Perform the equivalent of bash's touch on the file.

    Args:
        filename (str): File path to touch.
    '''

    open(filename, 'a').close()
    os.chmod(filename, stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


class TestData(object):

    def __init__(self, filename, dest_dir):
        self.filepath = None
        try:
            self.filepath = resource_filename(Requirement.parse("dammit"), 
                                              "dammit/tests/test-data/"     + filename)
        except ResolutionError:
            pass
        if not self.filepath or not os.path.isfile(self.filepath):
            self.filepath = os.path.join(os.path.dirname(__file__), 
                                          'test-data', filename)
        shutil.copy(self.filepath, dest_dir)
        self.filepath = os.path.join(dest_dir, filename)
    
    def __enter__(self):
        return self.filepath

    def __exit__(self, exc_type, exc_value, traceback):
        try:
            os.remove(self.filepath)
        except OSError:
            pass
        if exc_type:
            return False


class TemporaryFile(object):

    def __init__(self, directory):
        self.filepath = os.path.join(directory, str(hash(self)))

    def __enter__(self):
        return self.filepath

    def __exit__(self, exc_type, exc_value, traceback):
        try:
            os.remove(self.filepath)
        except OSError:
            pass
        if exc_type:
            return False


class Move(object):

    def __init__(self, target):
        print('Move to', target, file=sys.stderr)
        self.target = target
   
    def __enter__(self):
        self.cwd = os.getcwd()
        print('cwd:', self.cwd, file=sys.stderr)
        os.chdir(self.target)

    def __exit__(self, exc_type, exc_value, traceback):
        os.chdir(self.cwd)
        if exc_type:
            return False


class TemporaryDirectory(object):
    """Create and return a temporary directory.  This has the same
    behavior as mkdtemp but can be used as a context manager.  For
    example:

        with TemporaryDirectory() as tmpdir:
            ...

    Upon exiting the context, the directory and everything contained
    in it are removed.

    Note:
        Taken from http://stackoverflow.com/questions/19296146/tempfile-temporarydirectory-context-manager-in-python-2-7
    """

    def __init__(self, suffix="", prefix="tmp", dir=None):
        self._closed = False
        self.name = None # Handle mkdtemp raising an exception
        self.name = mkdtemp(suffix, prefix, dir)

    def __repr__(self):
        return "<{} {!r}>".format(self.__class__.__name__, self.name)

    def __enter__(self):
        return self.name

    def cleanup(self, _warn=False):
        if self.name and not self._closed:
            try:
                self._rmtree(self.name)
            except (TypeError, AttributeError) as ex:
                # Issue #10188: Emit a warning on stderr
                # if the directory could not be cleaned
                # up due to missing globals
                if "None" not in str(ex):
                    raise
                print >>_sys.stderr, "ERROR: {!r} while cleaning up {!r}".format(ex, self,)
                return
            self._closed = True
            if _warn:
                self._warn("Implicitly cleaning up {!r}".format(self), ResourceWarning)

    def __exit__(self, exc, value, tb):
        self.cleanup()
        if exc is not None:
            return False

    def __del__(self):
        # Issue a ResourceWarning if implicit cleanup needed
        self.cleanup(_warn=True)

    # XXX (ncoghlan): The following code attempts to make
    # this class tolerant of the module nulling out process
    # that happens during CPython interpreter shutdown
    # Alas, it doesn't actually manage it. See issue #10188
    _listdir = staticmethod(os.listdir)
    _path_join = staticmethod(os.path.join)
    _isdir = staticmethod(os.path.isdir)
    _islink = staticmethod(os.path.islink)
    _remove = staticmethod(os.remove)
    _rmdir = staticmethod(os.rmdir)
    _warn = _warnings.warn

    def _rmtree(self, path):
        # Essentially a stripped down version of shutil.rmtree.  We can't
        # use globals because they may be None'ed out at shutdown.
        for name in self._listdir(path):
            fullname = self._path_join(path, name)
            try:
                isdir = self._isdir(fullname) and not self._islink(fullname)
            except OSError:
                isdir = False
            if isdir:
                self._rmtree(fullname)
            else:
                try:
                    self._remove(fullname)
                except OSError:
                    pass
        try:
            self._rmdir(path)
        except OSError:
            pass


'''
These script running functions were taken from the khmer project:
https://github.com/dib-lab/khmer/blob/master/tests/khmer_tst_utils.py
'''

def scriptpath(scriptname='dammit'):
    "Return the path to the scripts, in both dev and install situations."

    # note - it doesn't matter what the scriptname is here, as long as
    # it's some dammit script present in this version of dammit.

    path = os.path.join(os.path.dirname(__file__), "../../bin")
    if os.path.exists(os.path.join(path, scriptname)):
        return path

    path = os.path.join(os.path.dirname(__file__), "../../../EGG-INFO/bin")
    if os.path.exists(os.path.join(path, scriptname)):
        return path

    for path in os.environ['PATH'].split(':'):
        if os.path.exists(os.path.join(path, scriptname)):
            return path


def _runscript(scriptname, sandbox=False):
    """
    Find & run a script with exec (i.e. not via os.system or subprocess).
    """

    import pkg_resources
    ns = {"__name__": "__main__"}
    ns['sys'] = globals()['sys']

    try:
        pkg_resources.get_distribution("dammit").run_script(scriptname, ns)
        return 0
    except pkg_resources.ResolutionError as err:
        if sandbox:
            path = os.path.join(os.path.dirname(__file__), "../sandbox")
        else:
            path = scriptpath()

        scriptfile = os.path.join(path, scriptname)
        if os.path.isfile(scriptfile):
            if os.path.isfile(scriptfile):
                exec(compile(open(scriptfile).read(), scriptfile, 'exec'), ns)
                return 0
        elif sandbox:
            raise nose.SkipTest("sandbox tests are only run in a repository.")

    return -1


def runscript(scriptname, args, in_directory=None,
              fail_ok=False, sandbox=False):
    """Run a Python script using exec().
    Run the given Python script, with the given args, in the given directory,
    using 'exec'.  Mimic proper shell functionality with argv, and capture
    stdout and stderr.
    When using :attr:`fail_ok`=False in tests, specify the expected error.
    """
    sysargs = [scriptname]
    sysargs.extend(args)
    cwd = os.getcwd()

    try:
        status = -1
        oldargs = sys.argv
        sys.argv = sysargs

        oldout, olderr = sys.stdout, sys.stderr
        sys.stdout = StringIO()
        sys.stdout.name = "StringIO"
        sys.stderr = StringIO()

        if in_directory:
            os.chdir(in_directory)
        else:
            in_directory = cwd

        try:
            print('running:', scriptname, 'in:', in_directory, file=oldout)
            print('arguments', sysargs, file=oldout)

            status = _runscript(scriptname, sandbox=sandbox)
        except nose.SkipTest:
            raise
        except SystemExit as e:
            status = e.code
        except:
            traceback.print_exc(file=sys.stderr)
            status = -1
    finally:
        sys.argv = oldargs
        out, err = sys.stdout.getvalue(), sys.stderr.getvalue()
        sys.stdout, sys.stderr = oldout, olderr

        os.chdir(cwd)

    if status != 0 and not fail_ok:
        #print(out)
        #print(err)
        assert False, (status, out, err)

    return status, out, err


def run_shell_cmd(cmd, fail_ok=False, in_directory=None):
    cwd = os.getcwd()
    if in_directory:
        os.chdir(in_directory)

    print('running: ', cmd)
    try:
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        (out, err) = p.communicate()

        out = out.decode('utf-8')
        err = err.decode('utf-8')

        if p.returncode != 0 and not fail_ok:
            print('out:', out)
            print('err:', err)
            raise AssertionError("exit code is non zero: %d" % p.returncode)

        return (p.returncode, out, err)
    finally:
        os.chdir(cwd)
