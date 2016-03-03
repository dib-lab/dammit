#!/usr/bin/env python
from __future__ import print_function

from flask.ext.zodb import ZODB
from flask.ext.zodb import Dict as zdict
import json
import pandas as pd

from .. import common
from .. import parsers


def get_database(app):
    db = getattr(g, '_database', None)
    if db is None:
        db = g._database = ZODB(app.config['ZODB_STORAGE'])
    return db
