# Referee's own exact helpers (independent of the author's code).
from fractions import Fraction
from itertools import product, combinations_with_replacement
import random

def monomials(dim, deg):
    """all exponent tuples alpha in N^dim with |alpha| <= deg"""
    out = []
    def rec(prefix, left, k):
        if k == dim:
            out.append(tuple(prefix)); return
        for a in range(left+1):
            rec(prefix+[a], left-a, k+1)
    rec([], deg, 0)
    return out

def rank_exact(rows):
    """exact rank over Q of an integer matrix (list of lists), fraction-free Bareiss-like elimination"""
    M = [list(map(int, r)) for r in rows]
    if not M: return 0
    nr, nc = len(M), len(M[0])
    r = 0
    for c in range(nc):
        piv = None
        for i in range(r, nr):
            if M[i][c] != 0:
                piv = i; break
        if piv is None: continue
        M[r], M[piv] = M[piv], M[r]
        pr = M[r]
        for i in range(r+1, nr):
            if M[i][c] != 0:
                f, g = M[i][c], pr[c]
                M[i] = [g*M[i][j] - f*pr[j] for j in range(nc)]
                # reduce by gcd to keep sizes small
                from math import gcd
                gg = 0
                for v in M[i]:
                    gg = gcd(gg, abs(v))
                if gg > 1:
                    M[i] = [v//gg for v in M[i]]
        r += 1
        if r == nr: break
    return r

def evalmat(S, deg):
    dim = len(S[0])
    mons = monomials(dim, deg)
    rows = []
    for y in S:
        row = []
        for a in mons:
            v = 1
            for yi, ai in zip(y, a):
                v *= yi**ai
            row.append(v)
        rows.append(row)
    return rows

def reg(S):
    """least d with every function on S a polynomial of degree <= d (exact)"""
    S = list(S)
    n = len(S)
    d = 0
    while True:
        if rank_exact(evalmat(S, d)) == n:
            return d
        d += 1

def width(S, u):
    vals = [sum(a*b for a, b in zip(y, u)) for y in S]
    return max(vals) - min(vals)

def omega_signed(S):
    dim = len(S[0])
    tot = Fraction(0)
    for eps in product((1, -1), repeat=dim):
        tot += width(S, eps)
    return tot / 2**dim

def MW(S):
    dim = len(S[0])
    tot = Fraction(0)
    for P in product((0, 1), repeat=dim):
        tot += width(S, P)
    return 2*tot / 2**dim
