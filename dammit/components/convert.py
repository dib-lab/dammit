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
from ..utils import touch

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


@convert.command(name='shmlast-to-gff3')
@click.argument('input_filename')
@click.argument('output_filename')
@click.option('--dbxref', default='shmlast-alignment',
              help='Database source to use in the GFF3 Dbxref'
                   ' field. You really should set this!')
def shmlast_to_gff3_cmd(input_filename, output_filename, dbxref):
    '''Given the CSV output from shmlast, convert it to GFF3 and
    save the results.
    \f

    Args:
        input_filename (str): The input CSV.
        output_filename (str): Destination for GFF3 output.
        database (str): Tag to use in the GFF3 `Dbxref` field.
    '''

    it = pd.read_csv(input_filename, chunksize=10000)
    writer = GFF3Writer(output_filename, converter=shmlast_to_gff3,
                        database=dbxref)

    try:
        for group in it:
            writer.write(group)
    except EmptyFile:
        touch(output_filename)


@convert.command(name='hmmscan-to-gff3')
@click.argument('input_filename')
@click.argument('output_filename')
@click.option('--dbxref', default='hmmscan-alignment',
              help='Database source to use in the GFF3 Dbxref'
                   ' field. You really should set this!')
def hmmscan_to_gff3_cmd(input_filename, output_filename, dbxref):
    '''Given HMMER output, convert it to GFF3 and
    save the results. The HMMER output can either be raw or
    dammit-converted zero-indexed, half-open interval CSV
    (as the good lord intended).
    \f

    Args:
        input_filename (str): The input CSV.
        output_filename (str): Destination for GFF3 output.
        database (str): Tag to use in the GFF3 `Dbxref` field.
    '''

    if input_filename.endswith('.csv') or input_filename.endswith('.tsv'):
        it = pd.read_csv(input_filename, chunksize=10000)
    else:
        it = HMMerParser(input_filename)

    writer = GFF3Writer(output_filename,
                        converter=hmmscan_to_gff3,
                        database=dbxref)
    try:
        for group in it:
            writer.write(group)
    except EmptyFile as e:
        touch(output_filename)


@convert.command(name='cmscan-to-gff3')
@click.argument('input_filename')
@click.argument('output_filename')
@click.option('--dbxref', default='hmmscan-alignment',
              help='Database source to use in the GFF3 Dbxref'
                   ' field. You really should set this!')
def cmscan_to_gff3_cmd(input_filename, output_filename, dbxref):
    '''Given raw input from Infernal's cmscan, convert it to GFF3 and
    save the results.
    \f

    Args:
        input_filename (str): The input CSV.
        output_filename (str): Destination for GFF3 output.
        database (str): Tag to use in the GFF3 `Dbxref` field.
    '''

    writer = GFF3Writer(output_filename, converter=cmscan_to_gff3,
                        database=dbxref)
    try:
        for group in InfernalParser(input_filename):
            writer.write(group)
    except EmptyFile as e:
        touch(output_filename)