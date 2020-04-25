# Copyright (C) 2015-2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the BSD license.  See the LICENSE file for details.

import argparse
import glob
import logging
import os
import sys
import subprocess
import yaml

from dammit import databases
from dammit import log
from dammit import utils
from dammit import ui
from dammit.meta import __version__, __authors__, __description__, __date__, __path__, get_config


class DammitApp(object):

    def handle_annotate(self):
        log.start_logging()
        print(ui.header('submodule: annotate', level=2))

        workflow_file = os.path.join(__path__, 'workflows', 'dammit.snakefile')
        config_file = os.path.join(__path__, 'config.yml')
        cmd = ["snakemake", "-s", workflow_file, "--configfile", config_file]

        if self.config_d['force'] is True:
            annot_targets, db_dir, out_dir = generate_targets(db=True, annot=True)
            utd_msg = '*All database tasks up-to-date.*'
            ood_msg = '*Some database tasks out-of-date; '\
                      'FORCE is True, ignoring!'
        #    uptodate, statuses = db_handler.print_statuses(uptodate_msg=utd_msg,
        #                                                   outofdate_msg=ood_msg)
        else:
            annot_targets, db_dir, out_dir = self.generate_targets(annot=True)
            #databases.check_or_fail(db_handler)

        #handle user-specified database dir and output dir properly
        # note if `--config` is last arg, it will try to add the workflow targets (targets) to config (and fail)
        config = ["--config", f"db_dir={db_dir}", f"dammit_dir={out_dir}", f"input_transcriptome={self.args.transcriptome}"]
        cmd.extend(config)

        helpful_args = ["-p", "--nolock", "--use-conda", "--rerun-incomplete", "-k", "--cores", f"{self.args.n_threads}"]
        cmd.extend(helpful_args)

        cmd.extend(annot_targets)

        print("Command: " + " ".join(cmd))
        try:
            subprocess.check_call(cmd)
        except subprocess.CalledProcessError as e:
            print(f'Error in snakemake invocation: {e}', file=sys.stderr)
            return e.returncode
