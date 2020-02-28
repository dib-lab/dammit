# Copyright (C) 2015-2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the BSD license.  See the LICENSE file for details.

import logging
import os
import sys


def default_database_dir(logger):
    '''Get the default database directory: checks the environment
    for a DAMMIT_DB_DIR variable, and if it is not found, returns
    the default location of `$HOME/.dammit/databases`.

    Args:
        logger (logging.logger): Logger to write to.
    Returns:
        str: Path to the database directory.
    '''

    try:
        directory = os.environ['DAMMIT_DB_DIR']
        logger.debug('found DAMMIT_DB_DIR env variable')
    except KeyError:
        logger.debug('no DAMMIT_DB_DIR or --database-dir, using'\
                     'default')
        directory = os.path.join(os.environ['HOME'], '.dammit', 'databases')
    return directory

