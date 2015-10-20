''' Optimized functions for remapping BLAST coordinates.

By default, BLAST outputs coordinates in 1-indexed fully open intervals.
Most programming languages and other programs use 0-index, half-open intervals.
This remapping can take a fair amount of time in pure python; with the help of
cython and numpy, it can be super fast.
'''

#!python
#cython: language_level=3, boundscheck=False

cimport numpy as np
import numpy as np


cdef np.ndarray[long] fix_coords_single(long sstart, long send, long qstart, long qend):
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
    cdef np.ndarray[long] res = np.empty(6, dtype=long)
    
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


cpdef np.ndarray[long, ndim=2] fix_blast_coords(np.ndarray[long] sstart, np.ndarray[long] send, 
                                                np.ndarray[long] qstart, np.ndarray[long] qend):
    '''Fix BLAST coordinates of many hits.

    Expects each column of attributes separately as numpy arrays.

    Args:
        sstart (ndarray): The subject start coordinates.
        send (ndarray): The subject end coordinates.
        qstart (ndarray): The query start coordinates.
        qend (ndarray): The query end coordinates.
    '''
    cdef long n = len(sstart)
    cdef long i = 0
    cdef np.ndarray[long, ndim=2] res = np.empty((n,6), dtype=long)
    for i in range(n):
        res[i,:] = fix_coords_single(sstart[i], send[i], qstart[i], qend[i])
    return res

