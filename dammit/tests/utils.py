#!/usr/bin/env python
from __future__ import print_function

import os as _os
import shutil
import stat
import sys
import warnings as _warnings
from pkg_resources import Requirement, resource_filename, ResolutionError

from tempfile import mkdtemp


def touch(filename):
    '''Perform the equivalent of bash's touch on the file.

    Args:
        filename (str): File path to touch.
    '''

    open(filename, 'a').close()
    _os.chmod(filename, stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


class TestData(object):

    def __init__(self, filename, dest_dir):
        self.filepath = None
        try:
            self.filepath = resource_filename(Requirement.parse("dammit"), 
                                              "dammit/tests/test-data/"     + filename)
        except ResolutionError:
            pass
        if not self.filepath or not _os.path.isfile(self.filepath):
            self.filepath = _os.path.join(_os.path.dirname(__file__), 
                                          'test-data', filename)
        shutil.copy(self.filepath, dest_dir)
        self.filepath = _os.path.join(dest_dir, filename)
    
    def __enter__(self):
        return self.filepath

    def __exit__(self, exc_type, exc_value, traceback):
        try:
            _os.remove(self.filepath)
        except OSError:
            pass
        if exc_type:
            return False


class TemporaryFile(object):

    def __init__(self, directory):
        self.filepath = _os.path.join(directory, str(hash(self)))

    def __enter__(self):
        return self.filepath

    def __exit__(self, exc_type, exc_value, traceback):
        try:
            _os.remove(self.filepath)
        except OSError:
            pass
        if exc_type:
            return False


class Move(object):

    def __init__(self, target):
        print('Move to', target, file=sys.stderr)
        self.target = target
   
    def __enter__(self):
        self.cwd = _os.getcwd()
        print('cwd:', self.cwd, file=sys.stderr)
        _os.chdir(self.target)

    def __exit__(self, exc_type, exc_value, traceback):
        _os.chdir(self.cwd)
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
    _listdir = staticmethod(_os.listdir)
    _path_join = staticmethod(_os.path.join)
    _isdir = staticmethod(_os.path.isdir)
    _islink = staticmethod(_os.path.islink)
    _remove = staticmethod(_os.remove)
    _rmdir = staticmethod(_os.rmdir)
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
