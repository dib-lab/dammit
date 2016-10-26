'''
Program metadata: the version, install path, description, and default config.
'''
import json
import os

__path__ = os.path.dirname(__file__)
__version__ = open(os.path.join(__path__, 'VERSION')).read().strip()
__authors__ = ['Camille Scott']
__description__ = 'a tool for easy de novo transcriptome annotation'
__date__ = 2016


def get_config():
    '''Parse the default JSON config files and return them as dictionaries.

    Returns:
        tuple: The config and databases dictionaries.
    '''
    with open(os.path.join(__path__, 'config.json')) as fp:
        config_d = json.load(fp)
    with open(os.path.join(__path__, 'databases.json'), 'r') as fp:
        databases_d = json.load(fp)
    return config_d, databases_d
