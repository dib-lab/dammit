#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2020
# File   : config.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 17.04.2020

from copy import deepcopy
import os
import yaml

from dammit.meta import __path__
from dammit.utils import read_yaml, update_nested_dict


def expand_user_config(config):
    for key, val in config.items():
        if key.endswith('dir'):
            config[key] = os.path.expandvars(val)


def merge_configs(default_user_config, internal_config,
                  external_config_file=None):

    config = deepcopy(default_user_config)
    if external_config_file:
        external_config = read_yaml(external_config_file)
        update_nested_dict(config, external_config)
        expand_user_config(config)
    update_nested_dict(config, internal_config)

    return config

#
# GLOBAL CONFIG VARIABLES
# 

DATABASES       = read_yaml(os.path.join(__path__, 'databases.yml'))
PIPELINES       = read_yaml(os.path.join(__path__, 'pipelines.yml'))
INTERNAL_CONFIG = read_yaml(os.path.join(__path__, 'internal.yml'))
USER_CONFIG     = read_yaml(os.path.join(__path__, 'config.yml'))

expand_user_config(USER_CONFIG)

def get_config_obj(external_config_file):

    class Obj:
        pass

    # Generate the config for *this specific run*
    CONFIG = Obj()
    CONFIG.databases = deepcopy(DATABASES)
    CONFIG.pipelines = deepcopy(PIPELINES)
    CONFIG.core = merge_configs(USER_CONFIG, INTERNAL_CONFIG, external_config_file)

    return CONFIG