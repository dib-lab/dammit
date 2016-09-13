#!/usr/bin/env python
from __future__ import print_function

from itertools import count 
import json
import logging
import os
import pprint
import re
from shutil import rmtree
import shutil
import sys

from doit.action import CmdAction
from doit.tools import run_once, create_folder, title_with_actions, LongRunning
from doit.task import clean_targets, dict_to_task

import pandas as pd
from khmer import HLLCounter, ReadParser
from shmlast.hits import BestHits
from shmlast.last import MafParser

from .utils import which, doit_task
from . import parsers
from . import gff


def clean_folder(target):
    try:
        rmtree(target)
    except OSError:
        pass


seq_ext = re.compile(r'(.fasta)|(.fa)|(.fastq)|(.fq)')
def strip_seq_extension(fn):
    return seq_ext.split(fn)[0]


@doit_task
def get_group_task(group_name, tasks):

    return {'name': group_name,
            'actions': None,
            'task_dep': [t.name for t in tasks]}

    

@doit_task
def get_download_task(url, target_fn):
    '''Creates a doit task to download the given URL.
    
    Args:
        url (str): URL to download.
        target_fn (str): Target for the download.
    Returns:
        dict: doit task.
    '''

    cmd = 'curl -o {target_fn} {url}'.format(**locals())
    name = 'download_gunzip:{0}'.format(os.path.basename(target_fn))

    return {'title': title_with_actions,
            'name': name,
            'actions': [LongRunning(cmd)],
            'targets': [target_fn],
            'clean': [clean_targets],
            'uptodate': [True]}



@doit_task
def get_download_and_gunzip_task(url, target_fn):
    '''Create a doit task which downloads and gunzips a file.

    Args:
        url (str): URL to download.
        target_fn (str): Target file for the download.
    Returns:
        dict: doit task.
    '''
    cmd = 'curl {url} | gunzip -c > {target_fn}'.format(**locals())

    name = 'download_and_gunzip:{0}'.format(os.path.basename(target_fn))

    return {'title': title_with_actions,
            'name': name,
            'actions': [LongRunning(cmd)],
            'targets': [target_fn],
            'clean': [clean_targets],
            'uptodate': [True]}



@doit_task
def get_download_and_untar_task(url, target_dir, label=None):
    '''Create a doit task to download a file and untar it in the
    given directory.

    Args:
        url (str): URL to download.
        target_dir (str: Directory to put the untarred folder in.
        label (str): Optional label to resolve doit name conflicts when putting
                     multiple results in the same folder.
    Returns:
        dict: doit task.
    '''

    if label is None:
        label = os.path.basename(url)

    cmd1 = 'mkdir -p {target_dir}; curl {url} | tar -xz -C {target_dir}'.format(**locals())
    name = 'download_and_untar:{0}-{1}'.format(os.path.basename(target_dir), label)
    done = os.path.join(target_dir, name) + '.done'
    cmd2 = 'touch {done}'.format(done=done)

    return {'name': name,
            'title': title_with_actions,
            'actions': [LongRunning(cmd1), cmd2],
            'targets': [done],
            'clean': [(clean_folder, [target_dir])],
            'uptodate': [True]}



@doit_task
def get_create_folder_task(folder):

    name = 'create_folder:{folder}'.format(**locals())

    return {'title': title_with_actions,
            'name': name,
            'actions': [(create_folder, [folder])],
            'targets': [folder],
            'uptodate': [run_once],
            'clean': [clean_targets] }


@doit_task
def get_sanitize_fasta_task(input_fn, output_fn):

    import re
    name = 'sanitize_fasta:{0}'.format(os.path.basename(input_fn))

    def fix():
        bad = r'["`/|=]+'
        with open(output_fn, 'wb') as fp:
            for record in ReadParser(input_fn):
                header = re.sub(bad, r'.', record.name)
                fp.write('>{0}\n{1}\n'.format(header, record.sequence))

    return {'name': name,
            'title': title_with_actions,
            'actions': [fix],
            'file_dep': [input_fn],
            'targets': [output_fn],
            'clean': [clean_targets]}


@doit_task
def get_rename_transcriptome_task(transcriptome_fn, output_fn, names_fn,
                                  transcript_basename, split_regex=None):

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
            'title': title_with_actions,
            'actions': [fix],
            'targets': [output_fn, names_fn],
            'file_dep': [transcriptome_fn],
            'clean': [clean_targets]}


@doit_task
def get_maf_best_hits_task(maf_fn, output_fn):

    hits_mgr = BestHits()

    def cmd():
        df = MafParser(maf_fn).read()
        df = hits_mgr.best_hits(df)
        df.to_csv(output_fn, index=False)

    name = 'maf_best_hits:{0}-{1}'.format(maf_fn, output_fn)

    return {'name': name,
            'title': title_with_actions,
            'actions': [cmd],
            'targets': [output_fn],
            'file_dep': [maf_fn],
            'clean': [clean_targets]}


@doit_task
def get_link_file_task(src, dst=''):
    ''' Soft-link file to the current directory
    '''
    cmd = 'ln -fs {src} {dst}'.format(src=src, dst=dst)
    return {'title': title_with_actions,
            'name': 'ln:' + os.path.basename(src) + ('-' + dst if dst else ''),
            'actions': [cmd],
            'file_dep': [src],
            'targets': [os.path.basename(src) if not dst else dst],
            'uptodate': [run_once],
            'clean': [clean_targets]}



@doit_task
def get_cat_task(file_list, target_fn):

    cmd = 'cat {files} > {t}'.format(files=' '.join(file_list), t=target_fn)

    return {'title': title_with_actions,
            'name': 'cat:' + os.path.basename(target_fn),
            'actions': [cmd],
            'file_dep': file_list,
            'targets': [target_fn],
            'clean': [clean_targets]}

   
@doit_task
def get_maf_gff3_task(input_filename, output_filename, database):

    name = 'maf-gff3:' + os.path.basename(output_filename)

    def cmd():
        if input_filename.endswith('.csv') or input_filename.endswith('.tsv'):
            it = pd.read_csv(input_filename, chunksize=10000)
        else:
            it = MafParser(input_filename)

        with open(output_filename, 'a') as fp:
            for group in it:
                gff_group = gff.maf_to_gff3_df(group, database=database)
                gff.write_gff3_df(gff_group, fp)

    return {'name': name,
            'title': title_with_actions,
            'actions': ['rm -f {0}'.format(output_filename),
                        cmd],
            'file_dep': [input_filename],
            'targets': [output_filename],
            'clean': [clean_targets]}


@doit_task
def get_crb_gff3_task(input_filename, output_filename, database):

    name = 'crbb-gff3:' + os.path.basename(output_filename)

    def cmd():
        with open(output_filename, 'a') as fp:
            for group in parsers.crb_to_df_iter(input_filename,
                                                remap=True):
                gff_group = gff.crb_to_gff3_df(group, database=database)
                gff.write_gff3_df(gff_group, fp)

    return {'name': name,
            'title': title_with_actions,
            'actions': ['rm -f {0}'.format(output_filename),
                        cmd],
            'file_dep': [input_filename],
            'targets': [output_filename],
            'clean': [clean_targets]}


@doit_task
def get_hmmscan_gff3_task(input_filename, output_filename, database):

    name = 'hmmscan-gff3:' + os.path.basename(output_filename)

    def cmd():
        with open(output_filename, 'a') as fp:
            for group in pd.read_csv(input_filename, chunksize=10000):
                gff_group = gff.hmmscan_to_gff3_df(group, database=database)
                gff.write_gff3_df(gff_group, fp)

    return {'name': name,
            'title': title_with_actions,
            'actions': ['rm -f {0}'.format(output_filename),
                        cmd],
            'file_dep': [input_filename],
            'targets': [output_filename],
            'clean': [clean_targets]}


@doit_task
def get_cmscan_gff3_task(input_filename, output_filename, database):

    name = 'cmscan-gff3:' + os.path.basename(output_filename)

    def cmd():
        with open(output_filename, 'a') as fp:
            for group in parsers.cmscan_to_df_iter(input_filename):
                gff_group = gff.cmscan_to_gff3_df(group, database=database)
                gff.write_gff3_df(gff_group, fp)

    return {'name': name,
            'title': title_with_actions,
            'actions': ['rm -f {0}'.format(output_filename),
                        cmd],
            'file_dep': [input_filename],
            'targets': [output_filename],
            'clean': [clean_targets]}


@doit_task
def get_gff3_merge_task(gff3_filenames, output_filename):

    name = 'gff3-merge:{0}'.format(os.path.basename(output_filename))

    merge_cmd = 'echo "{v}" > {out}; cat {f} | sed \'/^#/ d\''\
                ' | sort | sed \'/^$/d\' >> {out}'.format(v=gff.version_line(),
                                          f=' '.join(gff3_filenames),
                                          out=output_filename)
    return {'name': name,
            'title': title_with_actions,
            'actions': [merge_cmd],
            'file_dep': gff3_filenames,
            'targets': [output_filename],
            'clean': [clean_targets]}


@doit_task
def get_annotate_fasta_task(transcriptome_fn, gff3_fn, output_fn):

    name = 'fasta-annotate:{0}'.format(output_fn)

    def annotate_fasta():
        annotations = pd.concat([g for g in \
                                 parsers.parse_gff3(gff3_fn)])
        with open(output_fn, 'w') as fp:
            for n, record in enumerate(ReadParser(transcriptome_fn)):
                df = annotations.query('seqid == "{0}"'.format(record.name))
                annots = ['len={0}'.format(len(record.sequence))]
                #for seqid, sgroup in df.groupby('seqid'):
                #    for feature_type, fgroup in sgroup.groupby('feature_type'):
                for feature_type, fgroup in df.groupby('feature_type'):

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

                desc = '{0} {1}'.format(record.name, ' '.join(annots))

                fp.write('>{0}\n{1}\n'.format(desc.strip(), record.sequence))

    return {'name': name,
                'title': title_with_actions,
                'actions': [annotate_fasta],
                'file_dep': [transcriptome_fn, gff3_fn],
                'targets': [output_fn],
                'clean': [clean_targets]}


@doit_task
def get_transcriptome_stats_task(transcriptome, output_fn):

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
            'title': title_with_actions,
            'actions': [(cmd, [])],
            'file_dep': [transcriptome],
            'targets': [output_fn],
            'clean': [clean_targets]}
