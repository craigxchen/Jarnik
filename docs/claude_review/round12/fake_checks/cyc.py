"""Hilbert function of a T = 1 box invariant under the cyclic shift of the q = M-1 coordinates
(Paley-type profiles with base point inf = row 0, coordinates F_q), via isotypic decomposition mod p.

A mod-p rank is <= the rational rank; the Hilbert function mod p is therefore exact unless p divides some minor
(cross-checked with two primes).  The top degree computed mod p is >= rho_Q... (see write-up: we report both primes).
"""
import sys, time, itertools
import numpy as np
import flint
from sympy import isprime
from flib import *
from model import profile


def find_prime(q, start):
    p = start - (start % q) + 1
    while True:
        if p % q == 1 and isprime(p):
            return p
        p -= q


def root_of_unity(q, p):
    h = 2
    while True:
        w = pow(h, (p - 1) // q, p)
        if w != 1:
            return w
        h += 1


def shift(k, g, q):
    return tuple(k[(j - g) % q] for j in range(q))


def orbit_reps(items, q):
    seen = set()
    reps = []
    for k in items:
        if k in seen:
            continue
        orb = {shift(k, g, q) for g in range(q)}
        seen |= orb
        reps.append((min(orb), len(orb)))
    return reps


def check_invariant(pats, M):
    q = M - 1
    S = set(tuple(sorted(x - 1 for x in P)) for P in pats)
    S2 = set(tuple(sorted(((x + 1) % q) for x in P)) for P in S)
    return S == S2


def hilbert_cyc(A, q, maxdeg, p):
    preps = orbit_reps(A, q)
    mons = []
    for j in range(maxdeg + 1):
        mons += reduced_monomials(q, j)
    mreps = orbit_reps(mons, q)
    mreps.sort(key=lambda r: (sum(r[0]), r[0]))
    mdeg = [sum(a) for a, _ in mreps]
    w = root_of_unity(q, p)
    Pm = np.array([k for k, _ in preps], dtype=np.int64)
    Htot = {}
    for s in range(q):
        rows = [i for i, (k, sz) in enumerate(preps) if sz == q or s == 0]
        acc = np.zeros((len(rows), len(mreps)), dtype=object)
        for g in range(q):
            cg = pow(w, (-s * g) % q, p)
            # value of monomial a at shifted point shift(k,g)
            sh = np.array([shift(tuple(Pm[i]), g, q) for i in rows], dtype=np.int64)
            for jm, (a, _) in enumerate(mreps):
                v = np.ones(len(rows), dtype=np.int64)
                for x, e in enumerate(a):
                    if e:
                        v = v * sh[:, x] ** e
                acc[:, jm] = (acc[:, jm] + v * cg) % p
        Mm = flint.nmod_mat(len(rows), len(mreps), [int(x) % p for x in acc.ravel()], p)
        # rank profile by degree: rank of leading column blocks
        prev = 0
        ncols = 0
        for j in range(maxdeg + 1):
            ncols = sum(1 for dg in mdeg if dg <= j)
            sub = flint.nmod_mat(len(rows), ncols, [int(Mm[i, c]) for i in range(len(rows)) for c in range(ncols)], p)
            r = sub.rank()
            Htot[j] = Htot.get(j, 0) + (r - prev)
            prev = r
    return len(preps), len(mreps), [Htot[j] for j in range(maxdeg + 1)]


if __name__ == '__main__':
    name = sys.argv[1]; L = int(sys.argv[2]); maxdeg = int(sys.argv[3])
    M, pats = profile(name)
    q = M - 1
    print(name, 'shift-invariant:', check_invariant(pats, M), flush=True)
    A = box(M, pats, L)
    for start in (1 << 58, 1 << 59):
        p = find_prime(q, start)
        t = time.time()
        npo, nmo, H = hilbert_cyc(A, q, maxdeg, p)
        print(name, 'L=%d |A|=%d p=%d orbits=%d monorbits=%d H(0..%d)=%s sum=%d (%.1fs)' % (
            L, len(A), p, npo, nmo, maxdeg, H, sum(H), time.time() - t), flush=True)
