#!/usr/bin/env python
from __future__ import print_function

from .common import get_dammit_dir

import logging
import os
import sys
import textwrap

class LogFormatter(logging.Formatter):

    def __init__(self, width=90, padding=10):
        super(LogFormatter, self).__init__('')
        self.width = width
        self.padding = padding

    def do_wrap(self, msg, pad):
        wrapped = textwrap.wrap(msg, self.width - self.padding)
        return ['{0}{1}'.format(pad, ln) for ln in wrapped]

    def format(self, record):
        
        if record.levelno < 40:
            pad = ' ' * self.padding
            wrapped = self.do_wrap(record.msg, pad)
            res = '\n'.join(wrapped) + '\n'
        else:
            #pad = '{0}# '.format(' ' * (self.padding/2))
            #wrapped = self.do_wrap(record.msg, pad)
            extra = '[{0}:{1}]'.format(record.name, record.levelname)
            #extra = extra.rjust((self.width + len(extra)) - len(wrapped[0]))
            #wrapped[0] = '{0}{1}'.format(wrapped[0], extra)
            res = record.msg + extra

        return res


class DammitLogger(object):
    '''Set up logging for the dammit application. We insulate it
    in a class to let us choose to only activate it when the program itself
    is run, effectively keeping the tests and any use of the API from
    being noisy.
    '''

    def __init__(self):
        self.log_dir = os.path.join(get_dammit_dir(), 'log')
        try:
            os.makedirs(self.log_dir)
        except OSError:
            pass
        self.log_file = os.path.join(self.log_dir, 'dammit-all.log')

        self.config = { 'format': '%(asctime)s %(name)s:%(funcName)s:%(lineno)d '\
                                  '[%(levelname)s] \n%(message)s\n-----',
                        'datefmt': '%m-%d %H:%M:%S',
                        'filename':self.log_file,
                        'filemode': 'a' }

        # By default, only log errors (to the console)
        self.logger = logging.getLogger(__name__)
        noop = logging.NullHandler()
        self.logger.addHandler(noop)

    def run(self):
        logging.basicConfig(level=logging.DEBUG, **self.config)

        self.console = logging.StreamHandler(sys.stderr)
        self.console.setLevel(logging.INFO)
        self.formatter = LogFormatter()
        self.console.setFormatter(self.formatter)
        logging.getLogger('').addHandler(self.console)
        logging.getLogger('').debug('*** dammit! begin ***')


class LogReporter(object):
    """Log reporter. print results to the logger
    """
    # short description, used by the help system
    desc = 'console output'

    def __init__(self, logger):
        # save non-succesful result information (include task errors)
        self.failures = []
        self.runtime_errors = []
        self.logger = logger

    def write(self, text, level=logging.DEBUG):
        self.logger.log(level, text)

    def initialize(self, tasks):
        """called just after tasks have benn loaded before execution starts"""
        pass

    def get_status(self, task):
        """called when task is selected (check if up-to-date)"""
        pass

    def execute_task(self, task):
        """called when excution starts"""
        # ignore tasks that do not define actions
        # ignore private/hidden tasks (tasks that start with an underscore)
        if task.actions and (task.name[0] != '_'):
            self.write('[ ] %s' % task.name, logging.INFO)

    def add_failure(self, task, exception):
        """called when excution finishes with a failure"""
        self.failures.append({'task': task, 'exception':exception})

    def add_success(self, task):
        """called when excution finishes successfuly"""
        pass

    def skip_uptodate(self, task):
        """skipped up-to-date task"""
        if task.name[0] != '_':
            self.write("[x] %s" % task.name, level=logging.INFO)

    def skip_ignore(self, task):
        """skipped ignored task"""
        self.write("[!] %s" % task.name)

    def cleanup_error(self, exception):
        """error during cleanup"""
        self.logger.write(exeception.get_msg(), level=logging.ERROR)

    def runtime_error(self, msg):
        """error from doit (not from a task execution)"""
        # saved so they are displayed after task failures messages
        self.runtime_errors.append(msg)

    def teardown_task(self, task):
        """called when starts the execution of teardown action"""
        pass

    def complete_run(self):
        """called when finshed running all tasks"""
        # if test fails print output from failed task
        for result in self.failures:
            self.write('Some tasks failed!', level=logging.ERROR)
            msg = '%s - taskid:%s' % (result['exception'].get_name(),
                                        result['task'].name)
            self.write(msg, level=logging.ERROR)
            self.write(result['exception'].get_msg(), level=logging.ERROR)
            task = result['task']
            
            out = "".join([a.out for a in task.actions if a.out])
            self.write("%s" % out, level=logging.DEBUG)
            
            err = "".join([a.err for a in task.actions if a.err])
            self.write("%s" % err, level=logging.DEBUG)

        if self.runtime_errors:
            self.write("Execution aborted.", level=logging.ERROR)
            self.write(" ".join(self.runtime_errors), level=logging.ERROR)
