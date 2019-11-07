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

import click
from doit.tools import run_once, create_folder, LongRunning
from doit.task import clean_targets, dict_to_task
from khmer import HLLCounter, ReadParser
import pandas as pd

from dammit.cli import component


seq_ext = re.compile(r'(.fasta)|(.fa)|(.fastq)|(.fq)')
def strip_seq_extension(fn):
    return seq_ext.split(fn)[0]


@component.command()
@click.argument('transcriptome_fn')
@click.argument('output_fn')
@click.argument('names_fn')
@click.option('--basename', default='Transcript', show_default=True)
@click.option('--split-regex', default=None, type=str)
def rename_transcriptome(transcriptome_fn,
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
        for record in ReadParser(transcriptome_fn):
            header = header_func(record.name)
            fp.write('>{0}\n{1}\n'.format(header, record.sequence))
            names.append((record.name, header))

    pd.DataFrame(names, columns=['original', 'renamed']).to_csv(names_fn,
                                                            index=False)


@component.command()
@click.argument('transcriptome_fn')
@click.argument('output_fn')
@click.option('-K', default=25, type=int, 
              help='k-mer size for k-mer analysis.')
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
                                   'Offending transcript: \n>{0}\n{1}\nbad'\
                                   .format(contig.name, contig.sequence))
            if 'N' in sequence:
                sequence = sequence.replace('N', 'A')
                n_ambiguous += 1

            hll.consume_string(sequence)
            gc_len += contig.sequence.count('C')
            gc_len += contig.sequence.count('G')
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

    exp_kmers = (lens - (k+1)).sum()
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

