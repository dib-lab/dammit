#!/usr/bin/env python
import numpy as np
import csv
import sys
import pandas as pd


def fix_coords_single(sstart, send, qstart, qend):
    '''Fix BLAST coordinates for the parameters from a singe hit.

    Args:
        sstart (long): Subject start.
        send (long): Subject end.
        qstart (long): Query start.
        qend (long): Query end.
    Returns:
        ndarray: Numpy array with values [sstart, send, qstart, qend, sstrand, qstrand],
        where sstrand and qstrand are the subject and query strands as a 1 or -1.
    '''
    res = np.empty(6, dtype=long)
    
    if sstart < send:
        res[0] = sstart - 1
        res[1] = send
        res[4] = 1
    else:
        res[0] = send
        res[1] = sstart + 1
        res[4] = -1
    
    if qstart < qend:
        res[2] = qstart - 1
        res[3] = qend
        res[5] = 1
    else:
        res[2] = qend
        res[3] = qstart + 1
        res[5] = -1
    
    return res


def fix_blast_coords(sstart, send, qstart, qend):
    '''Fix BLAST coordinates of many hits.

    Expects each column of attributes separately as numpy arrays.

    Args:
        sstart (ndarray): The subject start coordinates.
        send (ndarray): The subject end coordinates.
        qstart (ndarray): The query start coordinates.
        qend (ndarray): The query end coordinates.
    '''
    n = len(sstart)
    i = 0
    res = np.empty((n,6), dtype=long)
    for i in range(n):
        res[i,:] = fix_coords_single(sstart[i], send[i], qstart[i], qend[i])
    return res


def remap_blast_coords_df(df):
    coords = fix_blast_coords(df.sstart.values, df.send.values, df.qstart.values, df.qend.values)
    df['sstart'] = coords[:,0]
    df['send'] = coords[:,1]
    df['qstart'] = coords[:,2]
    df['qend'] = coords[:,3]
    df['sstrand'] = coords[coords:,4]
    df['qstrand'] = coords[:,5]
