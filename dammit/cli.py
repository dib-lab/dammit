# Copyright (C) 2015-2019 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the BSD license.  See the LICENSE file for details.

import logging
import os
import sys
import yaml

import click
from click import get_current_context

from dammit.log import start_logging
from dammit.config import Config, parse_config


CONFIG = Config(*parse_config())


class ShortChoice(click.Choice):

    def get_metavar(self, param):
        return f"[{'|'.join(self.choices[:5])}|...]"


@click.group()
@click.pass_context
def component(ctx):
    logger = logging.getLogger('dammit.component')
    start_logging()
    ctx.obj = CONFIG


from .components import *



