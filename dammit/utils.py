#!/usr/bin/env python

import os

class Move(object):

    def __init__(self, target, create=False):
        print('Move to', target, file=sys.stderr)
        self.target = target
        self.create = create
   
    def __enter__(self):
        self.cwd = os.getcwd()
        print('cwd:', self.cwd, file=sys.stderr)
        if self.create:
            try:
                os.mkdir(self.target)
            except OSError:
                pass
        os.chdir(self.target)

    def __exit__(self, exc_type, exc_value, traceback):
        os.chdir(self.cwd)
        if exc_type:
            return False


def which(program):
    '''Checks whether the given program (or program path) is valid and
    executable.

    NOTE: Sometimes copypasta is okay! This function came from stackoverflow:

        http://stackoverflow.com/a/377028/5109965

    Args:
        program (str): Either a program name or full path to a program.

    Returns:
        Return the path to the executable or None if not found
    '''
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


