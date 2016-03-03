#!/usr/bin/env python
from __future__ import print_function

import json
import logging
import os

from flask import Flask, Blueprint, current_app, abort, jsonify
from flask import render_template, redirect, url_for
from jinja2 import TemplateNotFound
import pandas as pd

from .. import common
from .. import parsers
from . import genomed3plot as gd3
from .database import db

from . import static_folder, template_folder

views = Blueprint('transcript_pane', __name__,
                  static_folder=static_folder,
                  template_folder=template_folder)

@views.errorhandler(404)
def error_404_page(error):
    return render_template('404.html', msg=error.description), 404

@views.route('/tracks/<transcript>')
def tracks(transcript):
    try:
        transcript_info = db['transcripts'][transcript]
    except KeyError:
        abort(404, 'No transcript named {0} in the annotation database :('.format(transcript))
    else:
        tracks = gd3.create_tracks(transcript, transcript_info['annotations'])
        length = transcript_info['length']
        result = {'data': {'tracks': tracks, 'length': length}}
    return jsonify(result)

@views.route('/transcripts')
def transcript_list():
    return jsonify({'transcripts': db['annotations'].keys()})

@views.route('/transcript/<transcript>')
def transcript(transcript):
    print (transcript)
    return render_template('index.html',
                           name=current_app.config['DIRECTORY'],
                           transcript=transcript)
