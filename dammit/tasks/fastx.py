# Copyright (C) 2015-2018 Camille Scott
# All rights reserved.
#
# This software may be modified and distributed under the terms
# of the BSD license.  See the LICENSE file for details.

from itertools import count
import json
import os
import re

from doit.tools import run_once, create_folder, LongRunning
from doit.task import clean_targets, dict_to_task
from khmer import HLLCounter, ReadParser
import pandas as pd

from dammit.fileio.gff3 import GFF3Parser
from dammit.profile import profile_task
from dammit.utils import which, doit_task


seq_ext = re.compile(r'(.fasta)|(.fa)|(.fastq)|(.fq)')
def strip_seq_extension(fn):
    return seq_ext.split(fn)[0]

    
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

    import re
    name = os.path.basename(transcriptome_fn)

    if split_regex is None:
        counter = count()
        header_func = lambda name: '{0}_{1}'.format(transcript_basename, next(counter))
    else:
        def header_func(header):
            results = re.search(split_regex, header).groupdict()
            try:
                header = results['name']
            except KeyError as err:
                err.message = 'Header regex should have a name field!'
                raise
            return header

    def fix():
        names = []
        with open(output_fn, 'w') as fp:
            for record in ReadParser(transcriptome_fn):
                header = header_func(record.name)
                fp.write('>{0}\n{1}\n'.format(header, record.sequence))
                names.append((record.name, header))
        pd.DataFrame(names, columns=['original', 'renamed']).to_csv(names_fn,
                                                                    index=False)

    return {'name': name,
            'actions': [fix],
            'targets': [output_fn, names_fn],
            'file_dep': [transcriptome_fn],
            'clean': [clean_targets]}


@doit_task
@profile_task
def get_transcriptome_stats_task(transcriptome, output_fn):
    '''Create a doit task to run basic metrics on a transcriptome.

    Args:
        transcriptome (str): The input FASTA file.
        output_fn (str): File to store the results.

    Returns:
        dict: A doit task.
    '''

    import re

    name = 'transcriptome_stats:' + os.path.basename(transcriptome)
    K = 25
    DNA = re.compile(r'[ACGTN]*$', flags=re.IGNORECASE)

    def parse(fn):
        hll = HLLCounter(.01, K)
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

    def cmd():
        lens, uniq_kmers, gc_perc, n_amb = parse(transcriptome)

        exp_kmers = (lens - (K+1)).sum()
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
                 '25_mers': int(exp_kmers),
                 '25_mers_unique': uniq_kmers,
                 'n_ambiguous': n_amb,
                 'redundancy': redundancy,
                 'GCperc': float(gc_perc)}
        
        with open(output_fn, 'w') as fp:
            json.dump(stats, fp, indent=4)

        with open(output_fn, 'r') as fp:
            print(fp.read())

    return {'name': name,
            'actions': [(cmd, [])],
            'file_dep': [transcriptome],
            'targets': [output_fn],
            'clean': [clean_targets]}
