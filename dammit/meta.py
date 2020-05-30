# Copyright (C) 2015-2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the BSD license.  See the LICENSE file for details.

'''
Program metadata: the version, install path, description, and default config.
'''

import datetime
import json
import yaml
import os
import time

TIME_FMT = '%Y-%m-%d.%H-%M-%S'

__path__ = os.path.dirname(__file__)
__version__ = open(os.path.join(__path__, 'VERSION')).read().strip()
__authors__ = ['Camille Scott', "N. Tessa Pierce"]
__description__ = 'a tool for easy de novo transcriptome annotation'
__year__ = datetime.datetime.today().year
__time__ = time.strftime(TIME_FMT)