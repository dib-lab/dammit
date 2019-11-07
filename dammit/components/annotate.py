# Copyright (C) 2015-2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the BSD license.  See the LICENSE file for details.

import os
import sys

import click
from khmer import ReadParser

from ..cli import component
from ..fileio.gff3 import GFF3Parser


def generate_sequence_name(original_name, sequence, annotation_df):
    pass


def generate_sequence_summary(original_name, sequence, annotation_df):
    '''Given a FASTA sequence's original name, the sequence itself,
    and a DataFrame with its corresponding GFF3 annotations, generate
    a summary line of the annotations in key=value format.

    Args:
        original_name (str): Original name of the sequence.
        sequence (str): The sequence itself.
        annotation_df (DataFrame): DataFrame with GFF3 format annotations.

    Returns:
        str: The new summary header.
    '''

    annots = ['len={0}'.format(len(sequence))]
    for feature_type, fgroup in annotation_df.groupby('type'):

        if feature_type in ['translated_nucleotide_match',
                            'protein_hmm_match',
                            'RNA_sequence_secondary_structure']:

            collapsed = ','.join(['{}:{}-{}'.format(row.Name.split(':dammit')[0],
                                                     int(row.start),
                                                     int(row.end)) \
                            for _, row in fgroup.iterrows()])
            if feature_type == 'translated_nucleotide_match':
                key = 'homologies'
            elif feature_type == 'protein_hmm_match':
                key = 'hmm_matches'
            else:
                key = 'RNA_matches'
            annots.append('{0}={1}'.format(key, collapsed))

        elif feature_type in ['exon', 'CDS', 'gene',
                              'five_prime_UTR', 'three_prime_UTR',
                              'mRNA']:

            collapsed = ','.join(['{}-{}'.format(int(row.start),
                                                 int(row.end)) \
                            for _, row in fgroup.iterrows()])
            annots.append('{0}={1}'.format(feature_type, collapsed))

    desc = '{0} {1}'.format(original_name, ' '.join(annots))

    return desc


@component.command()
@click.argument('transcriptome_fn')
@click.argument('gff3_fn')
@click.argument('output_fn')
def annotate_fasta(transcriptome_fn, gff3_fn, output_fn):
    '''Annotate the headers in a FASTA file with its corresponding GFF3 file
    and place the resulting FASTA file in output_fn.
    \f

    Args:
        transcriptome_fn (str): Path to the FASTA file.
        gff3_fn (str): Path to the GFF3 annotations.
        output_fn (str): Path to store the resulting annotated FASTA.
    '''

    annotations = GFF3Parser(gff3_fn).read()
    with open(output_fn, 'w') as fp:
        for n, record in enumerate(ReadParser(transcriptome_fn)):
            df = annotations.query('seqid == "{0}"'.format(record.name))
            desc = generate_sequence_summary(record.name, record.sequence,
                                             df)
            fp.write('>{0}\n{1}\n'.format(desc.strip(), record.sequence))


