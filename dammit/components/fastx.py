#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2019
# File   : fastx.py
# License: BSD-3
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 24.10.2019


from itertools import count
# 
import json
import os
import re
import sys

import click
from khmer import HLLCounter, ReadParser
import pandas as pd

from .. import cloup
from ope.io.gff3 import GFF3Parser


seq_ext = re.compile(r'(.fasta)|(.fa)|(.fastq)|(.fq)')
def strip_seq_extension(fn):
    return seq_ext.split(fn)[0]


def rename_fasta(fasta_fn,
                 output_fn,
                 names_fn,
                 basename,
                 split_regex):
    '''Copy a FASTA file and rename the headers.
    \f

    Args:
        transcriptome_fn (str): The FASTA file.
        output_fn (str): Destination to copy to.
        names_fn (str): Destination to the store mapping from old to new names.
        basename (str): String to contruct new names from.
        split_regex (regex): Regex to split the input names with; must contain
            a `name` field.

    '''

    if split_regex is None:
        counter = count()
        header_func = lambda name: '{0}_{1}'.format(basename, next(counter))
    else:
        def header_func(header):
            results = re.search(split_regex, header).groupdict()
            try:
                header = results['name']
            except KeyError as err:
                err.message = 'Header regex should have a name field!'
                raise
            return header

    names = []
    with open(output_fn, 'w') as fp:
        for record in ReadParser(fasta_fn):
            header = header_func(record.name)
            fp.write('>{0}\n{1}\n'.format(header, record.sequence))
            names.append((record.name, header))

    pd.DataFrame(names, columns=['original', 'renamed']).to_csv(names_fn,
                                                                index=False)


@click.command('rename-fasta')
@click.argument('fasta_fn')
@click.argument('output_fn')
@click.argument('names_fn')
@click.option('--basename', default='Transcript', show_default=True)
@click.option('--split-regex', default=None, type=str)
def rename_fasta_cmd(fasta_fn,
                     output_fn,
                     names_fn,
                     basename,
                     split_regex):
    ''' Copy a FASTA file and rename the headers.
    '''

    allowed = r'[a-zA-Z0-9_\-:|\.]+'
    if not re.fullmatch(allowed, basename):
        print('ERROR: --basename must conform to {allowed}, please simplify it.', file=sys.stderr)
        sys.exit(1)

    if split_regex is not None and basename != 'Transcript':
        print('NOTE: --split-regex supersedes --basename', file=sys.stderr)

    rename_fasta(fasta_fn, output_fn, names_fn, basename, split_regex)


def transcriptome_stats(transcriptome_fn, output_fn, k):
    '''Run basic metrics on a transcriptome.
    \f

    Args:
        transcriptome (str): The input FASTA file.
        output_fn (str): File to store the results.

    Returns:
        dict: A doit task.
    '''

    DNA = re.compile(r'[ACGTN]*$', flags=re.IGNORECASE)

    def parse(fn):
        hll = HLLCounter(.01, k)
        lens = []
        names = []
        gc_len = 0
        n_ambiguous = 0
        for contig in ReadParser(fn):
            sequence = contig.sequence
            lens.append(len(sequence))
            names.append(contig.name)

            if DNA.match(sequence) is None:
                raise RuntimeError('non-ACGTN characters not supported. '\
                                   'Offending transcript: \n>{0}\n{1}'\
                                   .format(contig.name, contig.sequence))
            
            contig_n_ambiguous = 0
            if 'N' in sequence:
                contig_n_ambiguous += sequence.count('N')
                n_ambiguous += contig_n_ambiguous
                sequence = sequence.replace('N', 'A')

            hll.consume_string(sequence)
            gc_len += sequence.count('C')
            gc_len += sequence.count('G')

            # just assume gc content is .5 as a prior i suppose
            gc_len += contig_n_ambiguous // 2

        S = pd.Series(lens, index=names)
        try:
            S.sort_values()
        except AttributeError:
            S.sort()
        gc_perc = float(gc_len) / S.sum()
        return S, hll.estimate_cardinality(), gc_perc, n_ambiguous

    def calc_NX(lens, X):
        N = lens.sum()
        threshold = (float(X) / 100.0) * N

        NXlen = 0
        NXpos = 0
        cms = 0
        for n, l in enumerate(lens):
            cms += l
            if cms >= threshold:
                NXlen = l
                NXpos = n
                break
        return NXlen, NXpos

    lens, uniq_kmers, gc_perc, n_amb = parse(transcriptome_fn)

    exp_kmers = (lens - k + 1).sum()
    redundancy = float(exp_kmers - uniq_kmers) / exp_kmers
    if redundancy < 0:
        redundancy = 0.0

    N50len, N50pos = calc_NX(lens, 50)
    stats = {'N': len(lens),
             'sum': int(lens.sum()),
             'min': int(lens.min()),
             'max': int(lens.max()),
             'med': int(lens.median()),
             'mean': float(lens.mean()),
             'N50len': int(N50len),
             'N50pos': int(N50pos),
             'k': k,
             'k_mers': int(exp_kmers),
             'k_mers_unique': uniq_kmers,
             'n_ambiguous': n_amb,
             'redundancy': redundancy,
             'GCperc': float(gc_perc)}

    with open(output_fn, 'w') as fp:
        json.dump(stats, fp, indent=4)

    with open(output_fn, 'r') as fp:
        print(fp.read())


@click.command('transcriptome-stats')
@click.argument('transcriptome_fn')
@click.argument('output_fn')
@click.option('-K', default=25, type=int, 
              help='k-mer size for k-mer analysis.')
def transcriptome_stats_cmd(transcriptome_fn, output_fn, k):
    ''' Run basic metrics on a transcriptome
    '''

    transcriptome_stats(transcriptome_fn, output_fn, k)


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


@click.command('annotate-fasta')
@click.argument('transcriptome_fn')
@click.argument('gff3_fn')
@click.argument('output_fn')
def annotate_fasta_cmd(transcriptome_fn, gff3_fn, output_fn):
    '''Annotate a FASTA file from a GFF3 file.
    '''

    annotate_fasta(transcriptome_fn, gff3_fn, output_fn)