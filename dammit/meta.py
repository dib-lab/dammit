# Copyright (C) 2015-2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the BSD license.  See the LICENSE file for details.

'''
Program metadata: the version, install path, description, and default config.
'''
import json
import yaml
import os

__path__ = os.path.dirname(__file__)
__version__ = open(os.path.join(__path__, 'VERSION')).read().strip()
__authors__ = ['Camille Scott', "N. Tessa Pierce"]
__description__ = 'a tool for easy de novo transcriptome annotation'
__date__ = 2020


def get_config():
    '''Parse the default YAML or JSON config files and return them as dictionaries.

    Returns:
        tuple: The config and databases dictionaries.
    '''
    with open(os.path.join(__path__, 'config.yml')) as fp:
        try:
            config_d = yaml.safe_load(fp) #json.load(fp)
        except yaml.YAMLError as exc:
            print(exc)
    with open(os.path.join(__path__, 'databases.yml'), 'r') as fp:
        try:
            databases_d = yaml.safe_load(fp) #json.load(fp)
        except yaml.YAMLError as exc:
            print(exc)
    with open(os.path.join(__path__, 'pipelines.yml'), 'r') as fp:
        try:
            pipelines_d = yaml.safe_load(fp) #json.load(fp)
        except yaml.YAMLError as exc:
            print(exc)
    return config_d, databases_d, pipelines_d


def get_databases():
    '''Parse the default YAML or JSON config files and return them as dictionaries.

    Returns:
        tuple: The config and databases dictionaries.
    '''
    with open(os.path.join(__path__, 'databases.yml'), 'r') as fp:
        try:
            databases_d = yaml.safe_load(fp) #json.load(fp)
        except yaml.YAMLError as exc:
            print(exc)
    return databases_d

