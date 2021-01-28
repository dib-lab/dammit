# Copyright (C) 2015-2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the BSD license.  See the LICENSE file for details.

import click
import pandas as pd

from ope.io import EmptyFile

from ope.io.infernal import InfernalParser
from ope.io.hmmer import HMMerParser
from ope.io.gff3 import GFF3Writer
from ope.io.maf import MafParser

from ..convert import CMScan_to_GFF3, HMMScan_to_GFF3, MAF_to_GFF3, Shmlast_to_GFF3
from ..utils import touch


@click.command(name='busco-to-gff3')
@click.argument('busco_table_filename')
@click.argument('transcript_lengths_filename')
@click.argument('output_filename')
def busco_to_gff3_cmd(busco_table_filename, transcript_lengths_filename, output_filename):
    '''Convert a BUSCO full_table.tsv to GFF3.

    Because BUSCO does not output exact coordinates in transcriptome
    mode, the features are set to cover the entire target transcript.
    If the BUSCO outputs are updated to include the exact coordinates in the
    future, this will be changed.
    \f

    Args:
        busco_table_filename (str): The input full_table.tsv.
        transcript_lengths_filename (str): A CSV file contained the lengths of all transcripts.
        output_filename (str): Destination for GFF3 output.
    '''

    busco_df = pd.read_csv(busco_table_filename, delimiter='\t',
                           names=['BUSCO_id', 'Status', 'Sequence', 'Score', 'Length'],
                           error_bad_lines=False,
                           comment='#')
    lens_df  = pd.read_csv(transcript_lengths_filename,
                           names=['ID', 'Length'],
                           header=0)

    print(busco_df)
    print(lens_df.head())

    merged_df = pd.merge(busco_df, lens_df, left_on='Sequence', right_on='ID')
    print(merged_df)




@click.command(name='maf-to-gff3')
@click.argument('input_filename')
@click.argument('output_filename')
@click.option('--dbxref', default='MAF-alignment',
              help='Database source to use in the GFF3 Dbxref'
                   ' field. You really should set this!')
def maf_to_gff3_cmd(input_filename, output_filename, dbxref):
    '''Convert MAF to GFF3.
    
    Given either a raw MAF file or a CSV file with the proper MAF
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
                        converter=MAF_to_GFF3,
                        database=dbxref)
    try:
        for group in it:
            writer.write(group)
    except EmptyFile:
        touch(output_filename)


@click.command(name='shmlast-to-gff3')
@click.argument('input_filename')
@click.argument('output_filename')
@click.option('--dbxref', default='shmlast-alignment',
              help='Database source to use in the GFF3 Dbxref'
                   ' field. You really should set this!')
def shmlast_to_gff3_cmd(input_filename, output_filename, dbxref):
    '''Convert shmlast CSV output to GFF3.
    \f

    Args:
        input_filename (str): The input CSV.
        output_filename (str): Destination for GFF3 output.
        database (str): Tag to use in the GFF3 `Dbxref` field.
    '''

    it = pd.read_csv(input_filename, chunksize=10000)
    writer = GFF3Writer(output_filename, converter=Shmlast_to_GFF3,
                        database=dbxref)

    try:
        for group in it:
            writer.write(group)
    except EmptyFile:
        touch(output_filename)


@click.command(name='hmmscan-to-gff3')
@click.argument('input_filename')
@click.argument('output_filename')
@click.option('--dbxref', default='hmmscan-alignment',
              help='Database source to use in the GFF3 Dbxref'
                   ' field. You really should set this!')
def hmmscan_to_gff3_cmd(input_filename, output_filename, dbxref):
    '''Convert HMMER to GFF3.
    
    Given HMMER output, convert it to GFF3 and
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
                        converter=HMMScan_to_GFF3,
                        database=dbxref)
    try:
        for group in it:
            writer.write(group)
    except EmptyFile as e:
        touch(output_filename)


@click.command(name='cmscan-to-gff3')
@click.argument('input_filename')
@click.argument('output_filename')
@click.option('--dbxref', default='hmmscan-alignment',
              help='Database source to use in the GFF3 Dbxref'
                   ' field. You really should set this!')
def cmscan_to_gff3_cmd(input_filename, output_filename, dbxref):
    '''Convert Infernal's cmscan output to GFF3.
    \f

    Args:
        input_filename (str): The input CSV.
        output_filename (str): Destination for GFF3 output.
        database (str): Tag to use in the GFF3 `Dbxref` field.
    '''

    writer = GFF3Writer(output_filename, converter=CMScan_to_GFF3,
                        database=dbxref)
    try:
        for group in InfernalParser(input_filename):
            writer.write(group)
    except EmptyFile as e:
        touch(output_filename)
