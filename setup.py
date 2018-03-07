#!/usr/bin/env python3
# Copyright (C) 2015-2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the BSD license.  See the LICENSE file for details.

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

if sys.version_info < (3, 5):
    print >> sys.stderr, "ERROR: dammit! requires python 3.5 or greater"
    sys.exit()

version = open('dammit/VERSION').read().strip()

def main():
    setup(  name = 'dammit',
            version = version,
            description = 'dammit!',
            url = 'https://github.com/camillescott/dammit',
            author = 'Camille Scott',
            author_email = 'camille.scott.w@gmail.com',
            license = 'BSD',
            test_suite = 'pytest-runner',
            tests_require = ['pytest',
                             'codecov'],
            packages = find_packages(),
            scripts = glob('bin/*'),
            install_requires = ['setuptools>=0.6.35',
                                'pandas>=0.18.1',
                                'khmer>=2.0',
                                'doit>=0.29.0',
                                'ficus>=0.1',
                                'matplotlib>=1.0',
                                'numexpr>=2.3.1',
                                'shmlast>=1.2',
                                'filelock'],
            include_package_data = True,
            zip_safe = False,  )

if __name__ == "__main__":
    main()
