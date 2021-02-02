# Copyright (C) 2015-2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the BSD license.  See the LICENSE file for details.

import os
from io import StringIO
import traceback
import stat
import subprocess
import sys


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


def touch(filename):
    '''Perform the equivalent of bash's touch on the file.

    Args:
        filename (str): File path to touch.
    '''

    open(filename, 'a').close()
    os.chmod(filename, stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


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


def _runscript(scriptname):
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
        path = scriptpath()

        scriptfile = os.path.join(path, scriptname)
        if os.path.isfile(scriptfile):
            if os.path.isfile(scriptfile):
                exec(compile(open(scriptfile).read(), scriptfile, 'exec'), ns)
                return 0

    return -1


def runscript(scriptname, args, in_directory=None,
              fail_ok=False):
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

            status = _runscript(scriptname)
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
        print('STATUS:', status)
        print('STDOUT:', out)
        print('STDERR:', err)

    return status, out, err
    

def run(*args, **kwargs):
    return runscript('dammit', args, **kwargs)


def run_shell_cmd(cmd, fail_ok=False, in_directory=None):
    cwd = os.getcwd()
    if in_directory:
        os.chdir(in_directory)

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
