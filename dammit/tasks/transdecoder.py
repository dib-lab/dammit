# Copyright (C) 2015-2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the BSD license.  See the LICENSE file for details.

import os

from doit.action import CmdAction
from doit.tools import LongRunning
from doit.task import clean_targets

from dammit.tasks.utils import clean_folder, DependentTask, InstallationError
from dammit.profile import profile_task
from dammit.utils import which, doit_task


class TransDecoderLongOrfsTask(DependentTask):

    def deps(self):
        longorfs = which('TransDecoder.LongOrfs')
        if longorfs is None:
            raise InstallationError('TransDecoder.LongOrfs not found.')
        if self.logger:
            self.logger.debug('TransDecoder.LongOrfs:' + longorfs)
        return longorfs

    @doit_task
    @profile_task
    def task(self, input_filename, params=None):
        '''Get a task to run `Transdecoder.LongOrfs`.

        Args:
            input_filename (str): FASTA file to analyze.
            params (list): Extra parameters to pass to the executable.

        Returns:
            dict: A doit task.
        '''

        name = 'TransDecoder.LongOrfs:' + os.path.basename(input_filename)

        exc = self.deps()
        cmd = [exc, '-t', input_filename]
        if params is not None:
            cmd.extend(params)
        cmd = ' '.join(cmd)

        return {'name': name,
                'actions': [cmd],
                'file_dep': [input_filename],
                'targets': [input_filename + '.transdecoder_dir/longest_orfs.pep'],
                'clean': [(clean_folder, [input_filename + '.transdecoder_dir'])]}


class TransDecoderPredictTask(DependentTask):

    def deps(self):
        predict = which('TransDecoder.Predict')
        if predict is None:
            raise InstallationError('TransDecoder.Predict not found.')
        else:
            if self.logger:
                logger.debug('TransDecoder.Predict:' + predict)
            return predict

    @doit_task
    @profile_task
    def task(self, input_filename, pfam_filename=None, params=None):
        '''Get a task to run `TransDecoder.Predict`.

        Args:
            input_filename (str): The FASTA file to analyze.
            pfam_filename (str): If HMMER has been run against Pfam, pass this
                file name to `--retain_pfam_hits`.
            params (list): Extra parameters to pass to the executable.

        Returns:
            dict: A doit task.
        '''

        name = 'TransDecoder.Predict:' + os.path.basename(input_filename)

        exc = self.deps()
        cmd = [exc, '-t', input_filename]
        file_dep = [input_filename,
                    input_filename + '.transdecoder_dir/longest_orfs.pep']
        if pfam_filename is not None:
            cmd.extend(['--retain_pfam_hits', pfam_filename])
            file_dep.append(pfam_filename)
        if params is not None:
            cmd.extend(params)
        cmd = ' '.join(cmd)

        return {'name': name,
                'actions': [cmd],
                'file_dep': file_dep,
                'targets': [input_filename + '.transdecoder' + ext \
                            for ext in ['.bed', '.cds', '.pep', '.gff3']],
                'clean': [clean_targets,
                         (clean_folder, [input_filename + '.transdecoder_dir'])]}

