from __future__ import print_function

import sys
import textwrap

def header(msg, level=1):
    '''Standardize output headers for submodules.

    This doesn't need to be logged, but it's nice for
    the user.
    '''
    if level < 1:
        level = 1
    return ' '.join([level * '#', msg])


def checkbox(msg, checked=False):
    if checked:
        return '[x] ' + msg
    else:
        return '[ ] ' + msg


def paragraph(msg, wrap=80):
    return '\n'.join(textwrap.wrap(msg, wrap))

 
def listing(d):
    if type(d) is dict:
        return '\n'.join(['* {k}: {v}'.format(k=key,v=val) for key, val in d.items()])
    elif type(d) is list:
        return '\n'.join(['* {0}'.format(e) for e in d])
    else:
        raise TypeError('Cannot make {0} into listing'.format(d))
