#!/usr/bin/env python
import csv
import sys
import numpy as np
import pandas as pd

from blast import remap_blast_coords_df as remap_blast

'''
dammit uses 0-based, half-open intervals as its internal representation, as
spake by Saint Dijkstra (structured programs be upon Him) on the 397868400th
integer of Our Unix and observed by all pious and good followers of the Code.
Beware, fellow bioinformaticians! Heisenbugs and untyped Lambda Calculus be
upon ye who break from this most holy writ! Wretched followers of the 1,
demon peddlers of the fully closed and fully open intervals, sewers of discord!
Doom and damnation, dammit!
'''


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
        try:
            df[c] = df[c].astype(dtypes[c])
        except KeyError:
            pass


def blast_to_df_iter(fn, delimiter='\t', chunksize=10000, remap=True):
    '''Iterator of DataFrames of length chunksize parsed from an
    NCBI BLAST+ `-outfmt6` file.

    Native BLAST+ uses an interval of the form [start,end), start >= 1. This
    changes to [end,start) when on the negative strand, apparently soley
    to make other bioinformaticians suffer.

    We convert to proper 0-based, half-open intervals.

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

    GFF3 uses a 1-based, fully closed interval. Truly the devil's format.

    We convert to proper 0-based, half-open intervals.

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

        # Repent, repent!
        gtf_df.start = gtf_df.start - 1

        yield gtf_df


def crb_to_df_iter(fn, chunksize=10000, remap=True):
    '''Iterator of DataFrames of length chunksize parsed from
    the results from CRBB version crb-blast 0.6.6.

    crb-blast is given the same treatment as BLAST+, as that's what
    it uses under the hood.

    We convert to proper 0-based, half-open intervals.

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

    HMMER uses 1-based, fully open intervals. Another format of the devil.

    We convert to proper 0-based, half-open intervals.

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
        df.hmm_coord_from = df.hmm_coord_from - 1
        df.ali_coord_from = df.ali_coord_from - 1
        df.env_coord_from = df.env_coord_from -1
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

    1-based, fully open intervals. Truly Infernal.

    We convert to proper 0-based, half-open intervals.

    Args:
        fn (str): Path to the cmscan file.
        chunksize (int): Hits per iteration.
    Yields:
        DataFrame: Pandas DataFrame with the cmscan hits.
    '''

    def build_df(data):
        df = pd.DataFrame(data, columns=[k for k, _ in cmscan_cols])
        convert_dtypes(df, dict(cmscan_cols))
        df.mdl_from = df.mdl_from - 1
        df.seq_from = df.seq_from - 1
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
