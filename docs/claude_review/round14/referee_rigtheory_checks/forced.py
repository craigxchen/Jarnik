import sys, ctypes, math, numpy as np
from glib import *
lib = ctypes.CDLL('./fdfsc.so')
def canon(v):
    for x in v:
        if x: return x > 0
def forced_set(ps, K, qmax):
    r = len(ps); V = [v for v in VK(r, K) if canon(v)]
    Va = np.array(V, dtype=np.int32).ravel(); F = []; frees = {}
    for q in primerange(3, qmax+1):
        if q in ps: continue
        m = q-1 if q % 4 == 1 else q+1; n = m//4
        if n == 1: F.append(q); continue   # n=1: every v collides
        pv = nu(ps, q)
        par = np.array(pv if n % 2 == 0 else [-1]*r, dtype=np.int32)
        reps = [d for d in range(1, n) if n % d == 0 and (n % 2 == 1 or d % 2 == pv[0])]
        out = np.zeros(r, dtype=np.int32)
        f = lib.isfree(r, n, len(V), Va.ctypes.data_as(ctypes.POINTER(ctypes.c_int)), par.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
                       len(reps), np.array(reps, dtype=np.int32).ctypes.data_as(ctypes.POINTER(ctypes.c_int)), out.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))
        if f:
            c = [int(x) for x in out]
            assert all(sum(a*b for a, b in zip(v, c)) % n for v in V)
            if n % 2 == 0: assert all(c[j] % 2 == pv[j] for j in range(r))
            frees[q] = c
        else: F.append(q)
    return F, frees
if __name__ == "__main__":
    ps = [int(x) for x in sys.argv[1].split(",")]; K = int(sys.argv[2]); qmax = int(sys.argv[3])
    F, fr = forced_set(ps, K, qmax)
    print("K", K, "ps", ps, "qmax", qmax, "#F", len(F), "max", max(F)); print("F", F)
