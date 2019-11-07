# Copyright (C) 2015-2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the BSD license.  See the LICENSE file for details.

import os

from dammit.components.fastx import (rename_transcriptome,
                                     transcriptome_stats)
from dammit.profile import profile_task
from dammit.utils import which, doit_task


@doit_task
def get_rename_transcriptome_task(transcriptome_fn, output_fn, names_fn,
                                  transcript_basename, split_regex=None):
    '''Create a doit task to copy a FASTA file and rename the headers.

    Args:
        transcriptome_fn (str): The FASTA file.
        output_fn (str): Destination to copy to.
        names_fn (str): Destination to the store mapping from old to new names.
        transcript_basename (str): String to contruct new names from.
        split_regex (regex): Regex to split the input names with; must contain
            a `name` field.

    Returns:
        dict: A doit task.
    '''

    name = os.path.basename(transcriptome_fn)

    return {'name': name,
            'actions': [(rename_transcriptome, [transcriptome_fn,
                                                output_fn,
                                                names_fn,
                                                transcript_basename,
                                                split_regex])],
            'targets': [output_fn, names_fn],
            'file_dep': [transcriptome_fn],
            'clean': [clean_targets]}


@doit_task
@profile_task
def get_transcriptome_stats_task(transcriptome, output_fn):
    name = 'transcriptome_stats:' + os.path.basename(transcriptome)
    return {'name': name,
            'actions': [(transcriptome_stats, [transcriptome,
                                               output_fn])],
            'file_dep': [transcriptome],
            'targets': [output_fn],
            'clean': [clean_targets]}
