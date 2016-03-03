#!/usr/bin/env python

import sys, platform

# Automatically download setuptools if not available
try:
    from setuptools import *
except ImportError:
    from distribute_setup import use_setuptools
    use_setuptools()
finally:
    from setuptools import *

from glob import glob

if sys.version_info < (2, 4):
    print >> sys.stderr, "ERROR: dammit! requires python 2.4 or greater"
    sys.exit()

import numpy

cmdclass = {}

version = open('dammit/VERSION').read().strip()

def main():
    setup(  name = 'dammit',
            version = version,
            description = 'dammit!',
            url = 'https://github.com/camillescott/dammit',
            author = 'Camille Scott',
            author_email = 'camille.scott.w@gmail.com',
            license = 'BSD',
            test_suite = 'nose.collector',
            tests_require = ['nose'],
            packages = find_packages(),
            scripts = glob('bin/*'),
            install_requires = ['setuptools>=0.6.35',
                                'pandas>=0.17',
                                'khmer>=2.0',
                                'doit>=0.29.0',
                                'Sphinx>1.3.1',
                                'sphinx-rtd-theme>=0.1.9',
                                'nose==1.3.4',
                                'nose-capturestderr==1.0',
                                'ficus>=0.1',
                                'matplotlib>=1.0',
                                'scikit-learn>=0.16',
                                'numexpr>=2.3.1'],
            include_package_data = True,
            zip_safe = False,  )

if __name__ == "__main__":
    main()
