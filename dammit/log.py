# Copyright (C) 2015-2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the BSD license.  See the LICENSE file for details.

import logging
import os
import sys
import textwrap


def init_default_logger():
    log_dir = os.path.join(os.environ['HOME'], '.dammit', 'log')
    try:
        os.makedirs(log_dir)
    except OSError:
        pass
    log_file = os.path.join(log_dir, 'dammit-all.log')
    log_fmt = '%(asctime)s %(name)s:%(funcName)s:%(lineno)d '\
              '[%(levelname)s] \n%(message)s\n-----'
    date_fmt = '%m-%d %H:%M:%S'
    config = { 'format': log_fmt,
                    'datefmt': date_fmt,
                    'filename': log_file,
                    'filemode': 'a' }

    # By default, only log errors (to the console)
    logger = logging.getLogger(__name__)
    noop = logging.NullHandler()
    logger.addHandler(noop)

    def run(filename=None, test=False):
        if filename is None:
            filename = log_file
        if test is True:
            filename = os.path.join(log_dir, 'dammit-tests.log')
            print('Logger in testing mode:', filename)
        logging.basicConfig(level=logging.DEBUG, **config)

        run_handler = logging.FileHandler(filename)
        run_handler.setLevel(logging.DEBUG)
        #formatter = LogFormatter()
        run_handler.setFormatter(logging.Formatter(fmt=log_fmt, datefmt=date_fmt))
        logging.getLogger('').addHandler(run_handler)

        logging.getLogger('').debug('*** dammit BEGIN ***')

    return run


start_logging = init_default_logger()


