#!/usr/bin/env python
import csv
import sys
import numpy as np
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

gff_cols = [('seqid', str),
            ('source', str),
            ('feature_type', str),
            ('start', float),
            ('end', float),
            ('score', float),
            ('strand', str),
            ('frame', float),
            ('attributes', str)]

crb_cols = [('query', str),
            ('subject', str),
            ('id', str),
            ('aln_len', int),
            ('evalue', float),
            ('bitscore', float),
            ('qrange', str),
            ('srange', str),
            ('qlen', int),
            ('slen', int)]


def convert_dtypes(df, dtypes):
    '''Convert the columns of a DataFrame to the types specified
    in the given dictionary, inplace.

    Args:
        df (DataFrame): The DataFrame to convert.
        dtypes (dict): Dictionary mapping columns to types.
    '''

    for c in df.columns:
        df[c] = df[c].astype(dtypes[c])


def blast_to_df_iter(fn, delimiter='\t', chunksize=10000, remap=False):
    '''Iterator of DataFrames of length chunksize parsed from an
    NCBI BLAST+ `-outfmt6` file.

    Args:
        fn (str): The results file.
        chunksize (int): Hits per iteration.
    Yields:
        DataFrame: Pandas DataFrme with the BLAST+ hits.
    '''

    for group in pd.read_table(fn, header=None, skipinitialspace=True,
                                names=[k for k, _ in blast_cols],
                                delimiter=delimiter, chunksize=chunksize):
        convert_dtypes(group, dict(blast_cols))
        if remap:
            remap_blast(group)
        yield group


def parse_busco_full(fn):
    '''Parses a BUSCO full result table into a Pandas DataFrame.

    Args:
        fn (str): The results file.
    Returns:
        DataFrame: The results DataFrame.
    '''

    df = pd.read_table(fn)
    return df.rename(columns={'#BUSCO_group': 'BUSCO_group'})


def parse_busco_summary(fn):
    '''Parses a BUSCO summary file into a JSON compatible
    dictionary.

    Args:
        fn (str): The summary results file.
    Returns:
        dict: The BUSCO results.
    '''

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


def parse_busco_multiple(fn_list, dbs=['metazoa', 'vertebrata']):
    '''Parses multiple BUSCO results summaries into an appropriately
    index DataFrame.

    Args:
        fn_list (list): List of paths to results files.
        dbs (list): List of BUSCO database names.
    Returns:
        DataFrame: The formated DataFrame.
    '''

    data = []
    for fn in fn_list:
        data.append(parse_busco_summary(fn))

    df = pd.DataFrame(data)
    df['fn'] = [os.path.basename(fn)[14:-14].strip('.') for fn in fn_list]
    df['db'] = None
    for db in dbs:
        idx = df.fn.str.contains(db)
        df.loc[idx,'db'] = db
        df.loc[idx,'fn'] = df.loc[idx, 'fn'].apply(lambda fn: fn[:fn.find(db)].strip('. '))

    return df


def parse_gff3(fn, chunksize=10000):
    '''Iterator over DataFrames of length chunksize from a given
    GTF/GFF file.

    Args:
        fn (str): Path to the file.
        chunksize (int): Rows per iteration.
    Yields:
        DataFrame: Pandas DataFrame with the results.
    '''

    def attr_col_func(col):
        d = {}
        for item in col.strip(';').split(';'):
            key, _, val = item.strip().partition('=')
            d[key] = val.strip('')
        return d


    # Read everything into a DataFrame
    for group in  pd.read_table(fn, delimiter='\t', comment='#', 
                                names=[k for k,_ in gff_cols], na_values='.',
                                converters={'attributes': attr_col_func},
                                chunksize=chunksize, header=None,
                                dtype=dict(gff_cols)):

        # Generate a new DataFrame from the attributes dicts, and merge it in
        gtf_df = pd.merge(group,
                          pd.DataFrame(list(group.attributes)),
                          left_index=True, right_index=True)
        del gtf_df['attributes']
    
        # Switch from [start, end] to [start, end)
        gtf_df.end = gtf_df.end + 1
        #convert_dtypes(gtf_df, dict(gff_cols))

        yield gtf_df


def crb_to_df_iter(fn, chunksize=10000, remap=False):
    '''Iterator of DataFrames of length chunksize parsed from
    the results from CRBB version crb-blast 0.6.6.

    Args:
        fn (str): The results file.
        chunksize (int): Hits per iteration.
    Yields:
        DataFrame: Pandas DataFrame with the CRBB hits.
    '''

    for group in pd.read_table(fn, header=None, names=[k for k, _ in crb_cols],
                                delimiter='\t', chunksize=chunksize):

        convert_dtypes(group, dict(crb_cols))

        qrange = group.qrange.str.partition('..')
        group['qstart'] = qrange[0].astype(int)
        group['qend'] = qrange[2].astype(int)
        del group['qrange']
        srange = group.srange.str.partition('..')
        group['sstart'] = srange[0].astype(int)
        group['send'] = srange[2].astype(int)
        del group['srange']


        if remap:
            remap_blast(group)
        yield group


def parse_busco(fn):
    '''Parse a single BUSCO summary file to a dictionary.

    Args:
        fn (str): The results file.
    Returns:
        dict: The parsed results.
    '''

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
    ''' Given a list of BUSCO results from different databases, produce
    an appropriately multi-indexed DataFrame of the results.

    Args:
        fn_list (list): The BUSCO summary files.
        dbs (list): The BUSCO databases used for these runs.
    Returns:
        DataFrame: The BUSCO results.
    '''

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


def hmmscan_to_df_iter(fn, chunksize=10000):
    '''Iterator over DataFrames of length chunksize from a given
    hmmscan result file.

    Args:
        fn (str): Path to the hmmscan file.
        chunksize (int): Hits per iteration.
    Yields:
        DataFrame: Pandas DataFrame with the hmmscan hits.
    '''
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
    '''Iterator over DataFrames of length chunksize from a given
    cmscan result file.

    Args:
        fn (str): Path to the cmscan file.
        chunksize (int): Hits per iteration.
    Yields:
        DataFrame: Pandas DataFrame with the cmscan hits.
    '''

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
    '''Iterator yeilding DataFrames of length chunksize
    from a given TransDecoder GFF file.

    Args:
        fn (str): Path to the TransDecoder gff file.
        chunksize (int): Rows per iteration.
    Yields:
        DataFrame: Pandas DataFrame with the predict gene features.
    '''

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
    '''Iterator yielding DataFrames of length chunksize holding MAF alignments.

    An extra column is added for bitscore, using the equation described here:
        http://last.cbrc.jp/doc/last-evalues.html

    Args:
        fn (str): Path to the MAF alignment file.
        chunksize (int): Alignments to parse per iteration.
    Yields:
        DataFrame: Pandas DataFrame with the alignments.
    '''

    def fix_sname(name):
        new, _, _ = name.partition(',')
        return new

    def build_df(data, LAMBDA, K):
        df = pd.DataFrame(data)
        df['s_name'] = df['s_name'].apply(fix_sname)
        setattr(df, 'LAMBDA', LAMBDA)
        setattr(df, 'K', K)
        df['bitscore'] = (LAMBDA * df['score'] - np.log(K)) / np.log(2)
        return df

    data = []
    LAMBDA = None
    K = None
    with open(fn) as fp:
        while (True):
            try:
                line = fp.next().strip()
            except StopIteration:
                break
            if not line:
                continue
            if line.startswith('#'):
                if 'lambda' in line:
                    meta = line.strip(' #').split()
                    meta = {k:v for k, _, v in map(lambda x: x.partition('='), meta)}
                    LAMBDA = float(meta['lambda'])
                    K = float(meta['K'])
                else:
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
                    if LAMBDA is None:
                        raise Exception("old version of lastal; please update")
                    yield build_df(data, LAMBDA, K)
                    data = []

    if data:
        yield build_df(data, LAMBDA, K)
    
