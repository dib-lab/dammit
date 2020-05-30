#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2020
# File   : config.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 17.04.2020

import os
import yaml

from dammit.meta import __path__

DEFAULT_CONFIG_DIR = os.path.join(os.environ['HOME'], '.config', 'dammit')
DEFAULT_DATABASES_DIR = os.path.join(DEFAULT_CONFIG_DIR, 'databases')
DEFAULT_TEMP_DIR = os.path.join(DEFAULT_CONFIG_DIR, 'temp')
CONDA_ENV_TEMPDIR = 'envs'
WORKFLOW_CONFIG_TEMPDIR = 'configs'
TEMP_SUBDIRS = [CONDA_ENV_TEMPDIR, WORKFLOW_CONFIG_TEMPDIR]


def parse_config():
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


def create_tempdirs(root):
    for child in TEMP_SUBDIRS:
        os.makedirs(os.path.join(root, child), exist_ok=True)


class Config:
    def __init__(self, core, databases, pipelines, logger=None):
        self.core = core
        self.databases = databases
        self.pipelines = pipelines
        self.logger = logger


#######
#
# GLOBAL CONFIG VAR
#
#######

CONFIG = Config(*parse_config())

