# Copyright (C) 2015-2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the BSD license.  See the LICENSE file for details.

import os

from doit.tools import run_once, create_folder, LongRunning
from doit.task import clean_targets, dict_to_task
import pandas as pd
from shmlast import hits

from dammit.utils import which, doit_task, touch
from dammit.fileio import EmptyFile
from dammit.fileio.maf import MafParser
from dammit.fileio.infernal import InfernalParser
from dammit.fileio.hmmer import HMMerParser
from dammit.fileio.gff3 import (GFF3Writer, maf_to_gff3, shmlast_to_gff3,
                           hmmscan_to_gff3, cmscan_to_gff3)


@doit_task
def get_maf_best_hits_task(maf_fn, output_fn):
    '''Doit task to get the best hits from a lastal MAF file.

    Args:
        maf_fn (str): Path to the MAF file.
        output_fn (str): Path to store resulting CSV file.

    Returns:
        dict: A doit task.
    '''

    hits_mgr = hits.BestHits()

    def cmd():
        # can write out an empty file
        df = MafParser(maf_fn).read()
        df = hits_mgr.best_hits(df)
        df.to_csv(output_fn, index=False)

    name = 'maf_best_hits:{0}-{1}'.format(maf_fn, output_fn)

    return {'name': name,
            'actions': [cmd],
            'targets': [output_fn],
            'file_dep': [maf_fn],
            'clean': [clean_targets]}


@doit_task
def get_maf_gff3_task(input_filename, output_filename, database):
    '''Given either a raw MAF file or a CSV file with the proper MAF
    colums, convert it to GFF3 and save the results.

    Args:
        input_filename (str): The input MAF or CSV.
        output_filename (str): Destination for GFF3 output.
        database (str): Tag to use in the GFF3 `Dbxref` field.

    Returns:
        dict: A doit task.
    '''

    name = 'maf-gff3:' + os.path.basename(output_filename)

    def cmd():
        if input_filename.endswith('.csv') or input_filename.endswith('.tsv'):
            it = pd.read_csv(input_filename, chunksize=10000)
        else:
            it = MafParser(input_filename)
        writer = GFF3Writer(output_filename, converter=maf_to_gff3,
                            database=database)
        try:
            for group in it:
                writer.write(group)
        except EmptyFile:
            touch(output_filename)

    return {'name': name,
            'actions': ['rm -f {0}'.format(output_filename),
                        cmd],
            'file_dep': [input_filename],
            'targets': [output_filename],
            'clean': [clean_targets]}


@doit_task
def get_shmlast_gff3_task(input_filename, output_filename, database):
    '''Given the CSV output from shmlast, convert it to GFF3 and
    save the results.

    Args:
        input_filename (str): The input CSV.
        output_filename (str): Destination for GFF3 output.
        database (str): Tag to use in the GFF3 `Dbxref` field.

    Returns:
        dict: A doit task.
    '''
    
    name = 'shmlast-gff3:' + os.path.basename(output_filename)
    
    def cmd():
        it = pd.read_csv(input_filename, chunksize=10000)
        writer = GFF3Writer(output_filename, converter=shmlast_to_gff3,
                            database=database)

        try:
            for group in it:
                writer.write(group)
        except EmptyFile:
            touch(output_filename)

    return {'name': name,
            'actions': ['rm -f {0}'.format(output_filename),
                        cmd],
            'file_dep': [input_filename],
            'targets': [output_filename],
            'clean': [clean_targets]}


@doit_task
def get_hmmscan_gff3_task(input_filename, output_filename, database):
    '''Given HMMER output converted to CSV, convert it to GFF3 and
    save the results. CSV generated from the DataFrame(s) returned by
    the HMMerParser.

    Args:
        input_filename (str): The input CSV.
        output_filename (str): Destination for GFF3 output.
        database (str): Tag to use in the GFF3 `Dbxref` field.

    Returns:
        dict: A doit task.
    '''

    name = 'hmmscan-gff3:' + os.path.basename(output_filename)

    def cmd():
        writer = GFF3Writer(output_filename, converter=hmmscan_to_gff3,
                            database=database)
        try:
            for group in pd.read_csv(input_filename, chunksize=10000):
                writer.write(group)
        except EmptyFile as e:
            touch(output_filename)
            
    return {'name': name,
            'actions': ['rm -f {0}'.format(output_filename),
                        cmd],
            'file_dep': [input_filename],
            'targets': [output_filename],
            'clean': [clean_targets]}


@doit_task
def get_cmscan_gff3_task(input_filename, output_filename, database):
    '''Given raw input from Infernal's cmscan, convert it to GFF3 and
    save the results.

    Args:
        input_filename (str): The input CSV.
        output_filename (str): Destination for GFF3 output.
        database (str): Tag to use in the GFF3 `Dbxref` field.

    Returns:
        dict: A doit task.
    '''

    name = 'cmscan-gff3:' + os.path.basename(output_filename)

    def cmd():
        writer = GFF3Writer(output_filename, converter=cmscan_to_gff3,
                            database=database)
        try:
            for group in InfernalParser(input_filename):
                writer.write(group)
        except EmptyFile as e:
            touch(output_filename)

    return {'name': name,
            'actions': ['rm -f {0}'.format(output_filename),
                        cmd],
            'file_dep': [input_filename],
            'targets': [output_filename],
            'clean': [clean_targets]}


@doit_task
def get_gff3_merge_task(gff3_filenames, output_filename):
    '''Given a list of GFF3 files, merge them all together.

    Args:
        gff3_filenames (list): Paths to the GFF3 files.
        output_filename (str): Path to pipe the results.

    Returns:
        dict: A doit task.
    '''

    name = 'gff3-merge:{0}'.format(os.path.basename(output_filename))

    merge_cmd = 'echo "{v}" > {out}; cat {f} | sed \'/^#/ d\''\
                ' | sort | sed \'/^$/d\' >> {out}'.format(v=GFF3Writer.version_line,
                                          f=' '.join(gff3_filenames),
                                          out=output_filename)
    return {'name': name,
            'actions': [merge_cmd],
            'file_dep': gff3_filenames,
            'targets': [output_filename],
            'clean': [clean_targets]}

