#!/usr/bin/env python

import sys, platform

try:
    from Cython.Distutils import build_ext
except ImportError:
    use_cython = False
else:
    use_cython = True

# Automatically download setuptools if not available
from distribute_setup import use_setuptools
use_setuptools()

from setuptools import *
from glob import glob

if sys.version_info < (2, 4):
    print >> sys.stderr, "ERROR: dammit! requires python 2.4 or greater"
    sys.exit()

import numpy

cmdclass = {}

def main():
    setup(  name = 'dammit',
            version = '0.0.4-alpha',
            description = 'dammit!',
            url = 'https://github.com/camillescott/dammit',
            author = 'Camille Scott',
            author_email = 'camille.scott.w@gmail.com',
            license = 'BSD',
            test_suite = 'nose.collector',
            packages = ['dammit'],
            scripts = glob('bin/*'),
            ext_modules = get_extension_modules(),
            install_requires = ['khmer', 'nose>=1.0', 'doit', 'pandas'],
            include_package_data = True,
            zip_safe = False,
            cmdclass = cmdclass  )

def get_extension_modules():
    extensions = []
    
    if use_cython:
        ext = '.pyx'
        cmdclass.update({ 'build_ext': build_ext })
    else:
        ext = '.c'

    extensions.append( Extension( "dammit.cblast", [ "dammit/cblast.pyx" ],
                                  include_dirs=['.'] + [numpy.get_include()] ) )
    return extensions     

def monkey_patch_numpy():
    try:
        import numpy
        numpy.test = None
    except:
        pass
        
if __name__ == "__main__":
    monkey_patch_numpy()
    main()
