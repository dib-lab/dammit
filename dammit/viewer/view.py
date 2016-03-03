#!/usr/bin/env python
from __future__ import print_function

import json
import logging
import os

from flask import Flask, Blueprint
from flask import render_template, redirect, url_for
from jinja2 import TemplateNotFound
import pandas as pd

from .. import common
from .. import parsers
from . import genomed3plot as gd3

static_folder = os.path.join(common.rel_path, 'viewer', 'static')
template_folder = os.path.join(common.rel_path, 'viewer', 'templates')
transcript_pane = Blueprint('transcript_pane', __file__)

class ViewHandler(Flask):

    def __init__(self, annotate_handler, **kwargs):
        self.directory = annotate_handler.directory
        self.annotation = annotate_handler


        self.data_fn = os.path.join(self.directory,
                                    self.annotation.final_gff3_fn)
        self.data_df = pd.concat(parsers.parse_gff3(self.data_fn))
        self.sources = self.data_df.source.unique()

        super(ViewHandler, self).__init__(__name__.split('.')[0],
                                     static_folder=static_folder,
                                     template_folder=template_folder,
                                     **kwargs)

        self.route('/')(self.index)
        self.route('/transcript_list.json')(self.transcript_list)
        self.route('/die')(self.die)
        self.route('/tracks/<transcript>')(self.tracks)

    def tracks(self, transcript):
        tracks = gd3.create_tracks(transcript, self.data_df, self.sources)


    def transcript_list(self):
        return json.dumps(list(set(self.data_df['seqid'])))

    def index(self):
        return render_template('index.html',
                               name=self.annotation.transcriptome_fn)

    def die(self):
        raise RuntimeError

    def handle(self):
        common.print_header('Starting viewer server', level=2)
        self.run(host='0.0.0.0', port=5001, debug=True)
