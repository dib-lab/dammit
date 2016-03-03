#!/usr/bin/env python

import logging
import log
logger = log.DammitLogger()

from .hits import BestHits
from .crbl import CRBL

import fileio

from . import parsers
from . import gff
from . import blast
from . import tasks

from . import annotate
from . import databases
from . import dependencies
from . import common
from . import report

import os
rel_path = os.path.dirname(__file__)
__version__ = open(os.path.join(rel_path, 'VERSION')).read().strip()
