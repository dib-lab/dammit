''' Optimized functions for remapping BLAST coordinates.

By default, BLAST outputs coordinates in 1-indexed fully open intervals.
Most programming languages and other programs use 0-index, hafl-open intervals.
This remapping can take a fair amount of time in pure python; with the help of
cython and numpy, it can be super fast.
'''
cimport numpy as np
import numpy as np

# Returns: 
# [0: sstart
#  1: send
#  2: qstart
#  3: qend
#  4: sstrand
#  5: qstrand]
cdef np.ndarray[long] fix_coords_single(long sstart, long send, long qstart, long qend):
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
    cdef long n = len(sstart)
    cdef long i = 0
    cdef np.ndarray[long, ndim=2] res = np.empty((n,6), dtype=long)
    for i in range(n):
        res[i,:] = fix_coords_single(sstart[i], send[i], qstart[i], qend[i])
    return res

