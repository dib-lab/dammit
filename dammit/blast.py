#!/usr/bin/env python
#import pyximport
import numpy
#pyximport.install(setup_args={"include_dirs":numpy.get_include()}, reload_support=True)
import cblast
import csv
import sys
import pandas as pd

def best_hits(df, comp_col='evalue'):
    assert type(df) is pd.DataFrame
    df['ind'] = df.index
    df.sort(columns=['ind', comp_col], inplace=True)
    df.drop_duplicates(subset='ind', inplace=True)
    del df['ind']

def get_orthologies(A, B, idx, comp_col='evalue'):

    assert type(A) is pd.DataFrame
    assert type(B) is pd.DataFrame
    assert type(idx) is pd.Index

    best_hits(A)
    best_hits(B)

    X = pd.merge(A, B, how='inner', left_on='sseqid', right_index=True)
    X = pd.merge(pd.DataFrame(index=idx),
                      X[(X.index == X.sseqid_y)], left_index=True, right_index=True, how='left')
    del X['sseqid_x']
    del X['sseqid_y']

    return X

def remap_blast_coords_df(df):
    coords = cblast.fix_blast_coords(df.sstart.values, df.send.values, df.qstart.values, df.qend.values)
    df['sstart'] = coords[:,0]
    df['send'] = coords[:,1]
    df['qstart'] = coords[:,2]
    df['qend'] = coords[:,3]
    df['sstrand'] = coords[:,4]
    df['qstrand'] = coords[:,5]
