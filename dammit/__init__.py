#!/usr/bin/env python

import logging
import log
logger = log.DammitLogger()

from .hits import BestHits
from .crbl import CRBL

import parsers
import gff
import blast
import tasks

import annotate
import databases
import dependencies
import common
import report

import os
rel_path = os.path.dirname(__file__)
__version__ = open(os.path.join(rel_path, 'VERSION')).read().strip()
