"""
Helper functions for FM Index
---


2023 - Younginn Park
"""
# BWT related functions
from karksand import *

def bwtFromSa(t, sa=None):
    """ BWT construction from SA """
    if sa is None:
        sa = direct_kark_sort(t)

    bw = []
    dollarRow = None
    for si in sa:
        if si == 0:
            dollarRow = len(bw)
            bw.append('$')
        else:
            bw.append(t[si - 1])

    return ''.join(bw), dollarRow

def reverseBwt(fm):
    ''' Make T from BWT(T) '''
    rowi = 0 # start in first row
    t = '$' # start with rightmost character
    while fm.bwt[rowi] != '$':
        c = fm.bwt[rowi]
        t = c + t # prepend to answer
        # jump to row that starts with c of same rank
        rowi = fm.first[c] + fm.cps.rank(fm.bwt, c, rowi) - 1
    return t

def partialReverseBwt(fm, rowi, length):
    ''' Make a partial reverse of T from BWT(T) '''
    # rowi - row index of the last character of the reversed part
    t = ''
    for _ in range(length):
        c = fm.bwt[rowi]
        if c == "$": break # don't loop around T
        t = c + t
        rowi = fm.first[c] + fm.cps.rank(fm.bwt, c, rowi) - 1
    return t


# Alignment related functions
import numpy as np

def smithWaterman(x, y, s):
    ''' Calculate local alignment values of sequences x and y using
        dynamic programming.  Return maximal local alignment value. '''
    V = np.zeros((len(x)+1, len(y)+1), dtype=int)
    for i in range(1, len(x)+1):
        for j in range(1, len(y)+1):
            V[i, j] = max(V[i-1, j-1] + s(x[i-1], y[j-1]), # diagonal
                          V[i-1, j  ] + s(x[i-1], '-'),    # vertical
                          V[i  , j-1] + s('-',    y[j-1]), # horizontal
                          0)                               # empty
    return V


def traceback_with_indices(x, y, s, V=None):
    """ Trace back from given cell in local-alignment matrix V and find aligned indices """
    if not V:
        V = smithWaterman(x, y, s)
    
    # get i, j for maximal cell
    i, j = np.unravel_index(np.argmax(V), V.shape)
    start_i, start_j = i, j
    end_i, end_j = i, j
    while (i > 0 or j > 0) and V[i, j] != 0:
        diag, vert, horz = 0, 0, 0
        if i > 0 and j > 0:
            diag = V[i-1, j-1] + s(x[i-1], y[j-1])
        if i > 0:
            vert = V[i-1, j] + s(x[i-1], '-')
        if j > 0:
            horz = V[i, j-1] + s('-', y[j-1])
        if diag >= vert and diag >= horz:
            match = x[i-1] == y[j-1]
            i -= 1; j -= 1
        elif vert >= horz:
            i -= 1
        else:
            j -= 1
        # Update end indices
        end_i, end_j = i, j

    return (end_i, start_i - 1), (end_j, start_j - 1)

def editCost(xc, yc):
    return 1 if xc == yc else -1

def cost(xc, yc):
    ''' Cost function: 2 to match, -6 to gap, -4 to mismatch '''
    if xc == yc: return 2 # match
    if xc == '-' or yc == '-': return -6 # gap
    return -4
