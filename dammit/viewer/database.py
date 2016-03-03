#!/usr/bin/env python
from __future__ import print_function

from flask import g, current_app
from flask.ext.zodb import ZODB

db = ZODB()
