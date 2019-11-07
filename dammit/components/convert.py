# Copyright (C) 2015-2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the BSD license.  See the LICENSE file for details.

import click
import pandas as pd

from ..cli import component
from ..fileio import EmptyFile

from ..fileio.infernal import InfernalParser
from ..fileio.hmmer import HMMerParser
from ..fileio.gff3 import (GFF3Writer, maf_to_gff3, shmlast_to_gff3,
                           hmmscan_to_gff3, cmscan_to_gff3)
from ..fileio.maf import MafParser

@component.group()
def convert():
    pass


@convert.command(name='maf-to-gff3')
@click.argument('input_filename')
@click.argument('output_filename')
@click.option('--dbxref', default='MAF-alignment',
              help='Database source to use in the GFF3 Dbxref'
                   ' field. You really should set this!')
def maf_to_gff3_cmd(input_filename, output_filename, dbxref):
    '''Given either a raw MAF file or a CSV file with the proper MAF
    colums, convert it to GFF3 and save the results.
    \f

    Args:
        input_filename (str): The input MAF or CSV.
        output_filename (str): Destination for GFF3 output.
        dbxref (str): Tag to use in the GFF3 `Dbxref` field.
    '''

    if input_filename.endswith('.csv') or input_filename.endswith('.tsv'):
        it = pd.read_csv(input_filename, chunksize=10000)
    else:
        it = MafParser(input_filename)
    writer = GFF3Writer(filename=output_filename,
                        converter=maf_to_gff3,
                        database=dbxref)
    try:
        for group in it:
            writer.write(group)
    except EmptyFile:
        touch(output_filename)



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

