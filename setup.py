#!/usr/bin/env python

from setuptools import setup
from os.path import join as path_join
from os import listdir as os_listdir


scripts = [path_join('bin', script) for script in os_listdir('bin')]

setup(name='dammit',
        version='0.0.1-alpha',
        description='dammit!',
        url='https://github.com/camillescott/dammit',
        author='Camille Scott',
        author_email='camille.scott.w@gmail.com',
        license='BSD',
        test_suite='nose.collector',
        packages=['dammit'],
        scripts=scripts,
        install_requires=['khmer', 'screed', 'nose>=1.0', 'doit', 'pandas'],
        include_package_data=True,
        zip_safe=False)
