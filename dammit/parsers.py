#!/usr/bin/env python
import csv
import sys
import pandas as pd

from blast import remap_blast_coords_df as remap_blast

blast_cols = [('qseqid', str), 
              ('sseqid', str), 
              ('pident', float), 
              ('length', int), 
              ('mismatch', int), 
              ('gapopen', int),
              ('qstart', int), 
              ('qend', int), 
              ('sstart', int), 
              ('send', int), 
              ('evalue', float), 
              ('bitscore', float)]



hmmscan_cols = [('target_name', str), 
                ('target_accession', str), 
                ('tlen', int),
                ('query_name', str), 
                ('query_accession', str), 
                ('query_len', int),
                ('full_evalue', float), 
                ('full_score', float), 
                ('full_bias', float), 
                ('domain_num', int), 
                ('domain_total', int), 
                ('domain_c_evalue', float), 
                ('domain_i_evalue', float),
                ('domain_score', float), 
                ('domain_bias', float), 
                ('hmm_coord_from', int), 
                ('hmm_coord_to', int), 
                ('ali_coord_from', int), 
                ('ali_coord_to', int), 
                ('env_coord_from', int), 
                ('env_coord_to', int), 
                ('accuracy', float), 
                ('description', str)]

cmscan_cols = [('target_name', str),
               ('target_accession', str),
               ('query_name', str),
               ('query_accession', str),
               ('mdl', str),
               ('mdl_from', int),
               ('mdl_to', int),
               ('seq_from', int),
               ('seq_to', int),
               ('strand', str), 
               ('trunc', str), 
               ('pass', str),
               ('gc', float),
               ('bias', float),
               ('score', float),
               ('e_value', float),
               ('inc', str),
               ('description', str)]

gff3_transdecoder_cols = [('seqid', str), 
                          ('feature_type', str), 
                          ('start', int), 
                          ('end', int), 
                          ('strand', str)]

def convert_dtypes(df, dtypes):
    for c in df.columns:
        df[c] = df[c].astype(dtypes[c])

def blast_to_df_iter(fn, delimiter='\t', chunksize=10000, remap=False):

    for group in  pd.read_table(fn, header=None, skipinitialspace=True,
                                names=[k for k, _ in blast_cols],
                                delimiter=delimiter):
        convert_dtypes(group, dict(blast_cols))
        if remap:
            remap_blast(group)
        yield group

def parse_busco(fn):
    res = {}
    with open(fn) as fp:
        for ln in fp:
            if ln.strip().startswith('C:'):
                tokens = ln.split(',')
                for token in tokens:
                    key, _, val = token.partition(':')
                    key = key.strip()
                    val = val.strip().strip('%')
                    if key == 'C':
                        valc, _, vald = val.partition('%')
                        valc = valc.strip()
                        vald = vald.strip('D:][%')
                        res['C(%)'] = valc
                        res['D(%)'] = vald
                    else:
                        if key != 'n':
                           key += '(%)'
                        res[key] = val.strip().strip('%')
    return res


def busco_to_df(fn_list, dbs=['metazoa', 'vertebrata']):

    data = []
    for fn in fn_list:
        data.append(parse_busco(fn))

    df = pd.DataFrame(data)
    df['fn'] = [os.path.basename(fn)[14:-14].strip('.') for fn in fn_list]
    df['db'] = None
    for db in dbs:
        idx = df.fn.str.contains(db)
        df.loc[idx,'db'] = db
        df.loc[idx,'fn'] = df.loc[idx, 'fn'].apply(lambda fn: fn[:fn.find(db)].strip('. '))
    return df


def gtf_to_df(filename):

    # Converter func for the nonstandard attributes column
    def attr_col_func(col):
        d = {}
        for item in col.strip(';').split(';'):
            pair = item.strip().split(' ')
            d[pair[0]] = pair[1].strip('"')
        return d

    def strand_func(col):
        if col =='+':
            return 1
        else:
            return -1

    names=['contig_id', 'source', 'feature', 'start', 'end',
           'score', 'strand', 'frame', 'attributes']

    # Read everything into a DataFrame
    gtf_df = pd.read_table(filename, delimiter='\t', comment='#',
                           header=False, names=names,
                           converters={'attributes': attr_col_func, 'strand': strand_func})
    
    # Generate a new DataFrame from the attributes dicts, and merge it in
    gtf_df = pd.merge(gtf_df,
                      pd.DataFrame(list(gtf_df.attributes)),
                      left_index=True, right_index=True)
    del gtf_df['attributes']
    
    # Switch from [start, end] to [start, end)
    gtf_df.end = gtf_df.end + 1

    return gtf_df


def hmmscan_to_df_iter(fn, chunksize=10000):
    
    def split_query(item):
        q, _, _ = item.rpartition('|')
        return q

    def build_df(data):
        df = pd.DataFrame(data, columns=[k for k, _ in hmmscan_cols])
        convert_dtypes(df, dict(hmmscan_cols))
        df['full_query_name'] = df.query_name
        df['query_name'] = df.query_name.apply(split_query)
        return df

    data = []
    with open(fn) as fp:
        for n, ln in enumerate(fp):
            if not ln or ln.startswith('#'):
                continue

            tokens = ln.split()
            data.append(tokens[:len(hmmscan_cols)-1] + \
                        [' '.join(tokens[len(hmmscan_cols)-1:])])
            if len(data) >= chunksize:
                yield build_df(data)
                data = []

    if data:
        yield build_df(data)

def cmscan_to_df_iter(fn, chunksize=10000):
    
    def build_df(data):
        df = pd.DataFrame(data, columns=[k for k, _ in cmscan_cols])
        convert_dtypes(df, dict(cmscan_cols))
        return df

    data = []
    with open(fn) as fp:
        for ln in fp:
            ln = ln.strip()
            if not ln or ln.startswith('#'):
                continue
            tokens = ln.split()
            data.append(tokens[:len(cmscan_cols)-1] + \
                        [' '.join(tokens[len(cmscan_cols)-1:])])

            if len(data) >= chunksize:
                yield build_df(data)
                data = []

    if data:
        yield build_df(data)

def gff3_transdecoder_to_df_iter(fn, chunksize=10000):

    def build_df(data):
        df = pd.DataFrame(data, columns=[k for k, _ in gff3_transdecoder_cols])
        convert_dtype(df, dict(gff3_transdecoder_cols))
        return df
    
    data = []
    with open(fn) as fp:
        for ln in fp:
            if ln == '\n':
                continue
            tokens = ln.split('\t')
            try:
                data.append([tokens[0]] + tokens[2:5] + [tokens[6]])
            except IndexError as e:
                print e
                print tokens
                break
            if len(data) >= chunksize:
                yield build_df(data)
                data = []

    if data:
        yield build_df(data)

def maf_to_df_iter(fn, chunksize=10000):
    
    def build_df(data):
        return pd.DataFrame(data)

    data = []
    with open(fn) as fp:
        while (True):
            try:
                line = fp.next().strip()
            except StopIteration:
                break
            if not line or line.startswith('#'):
                continue
            if line.startswith('a'):
                cur_aln = {}

                # Alignment info
                tokens = line.split()
                for token in tokens[1:]:
                    key, _, val = token.strip().partition('=')
                    cur_aln[key] = float(val)
                
                # First sequence info
                line = fp.next()
                tokens = line.split()
                cur_aln['s_name'] = tokens[1]
                cur_aln['s_start'] = int(tokens[2])
                cur_aln['s_aln_len'] = int(tokens[3])
                cur_aln['s_strand'] = tokens[4]
                cur_aln['s_len'] = int(tokens[5])
                
                # First sequence info
                line = fp.next()
                tokens = line.split()
                cur_aln['q_name'] = tokens[1]
                cur_aln['q_start'] = int(tokens[2])
                cur_aln['q_aln_len'] = int(tokens[3])
                cur_aln['q_strand'] = tokens[4]
                cur_aln['q_len'] = int(tokens[5])

                data.append(cur_aln)
                if len(data) >= chunksize:
                    yield build_df(data)
                    data = []

    if data:
        yield build_df(data)
    
