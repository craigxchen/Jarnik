"""Hilbert function of T = 1 slab-boxes with a cyclic symmetry, via isotypic decomposition mod p.

Coordinates indexed by Z_q (the non-base points), the shift x -> x+1 preserves A.  Columns: translates
i + D of a base pattern D (Paley: D = {0} u QR_q), slab |k(i + D)| <= L.  Optional extra columns.

For each character s of Z_q, the s-isotypic part of the function space on A is determined by values on orbit
representatives; the projections P_s(mu) of monomial-orbit representatives, in degree order, give the rank
profile, hence H_A(i) = sum_s H_s(i).  Mod-p ranks (p = a prime = 1 mod q); a mod-p rank is <= the rational rank,
so mod-p Hilbert values can only be smaller in low degrees and the TOP DEGREE rho computed mod p satisfies
rho_p >= rho_Q? -- no: see note in the write-up; we report the mod-p values as evidence and cross-check with two
primes.
"""
import itertools, sys, time
import numpy as np
import flint
from sympy import nextprime, isprime


def find_prime(q, start):
    p = start - (start % q) + 1
    while True:
        if p % q == 1 and isprime(p):
            return p
        p -= q


def root_of_unity(q, p):
    import random
    random.seed(1)
    while True:
        h = random.randrange(2, p - 1)
        w = pow(h, (p - 1) // q, p)
        if w != 1:
            return w


def points(q, D, L, extra=()):
    pts = []
    cols = [[(i + x) % q for x in D] for i in range(q)] + [list(S) for S, _ in extra]
    Ls = [L] * q + [Lx for _, Lx in extra]
    for k in itertools.product((-1, 0, 1), repeat=q):
        if all(abs(sum(k[x] for x in S)) <= Lc for S, Lc in zip(cols, Ls)):
            pts.append(k)
    return pts


def shift(k, g, q):
    # (T^g k)_j = k_{j-g}
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


def masks(vecs, q):
    """for ternary point vectors: nz mask, neg mask (int64 arrays)."""
    V = np.array(vecs, dtype=np.int64)
    pw = (1 << np.arange(q, dtype=np.int64))
    nz = ((V != 0).astype(np.int64) * pw).sum(axis=1)
    ng = ((V < 0).astype(np.int64) * pw).sum(axis=1)
    return nz, ng


def popcount_parity(x):
    x = x.copy()
    par = np.zeros_like(x)
    while True:
        nzm = x != 0
        if not nzm.any():
            break
        par ^= (x & 1)
        x >>= 1
    return par


def hilbert_cyclic(q, D, L, p, extra=(), verbose=True):
    t0 = time.time()
    pts = points(q, D, L, extra)
    N = len(pts)
    preps = orbit_reps(pts, q)
    mons = list(itertools.product((0, 1, 2), repeat=q))
    mreps = orbit_reps(mons, q)
    mreps.sort(key=lambda r: (sum(r[0]), r[0]))
    mdeg = np.array([sum(a) for a, _ in mreps])
    A1 = np.array([[1 if x == 1 else 0 for x in a] for a, _ in mreps], dtype=np.int64)
    A2 = np.array([[1 if x == 2 else 0 for x in a] for a, _ in mreps], dtype=np.int64)
    pw = (1 << np.arange(q, dtype=np.int64))
    a1 = (A1 * pw).sum(axis=1)
    asupp = ((A1 + A2) * pw).sum(axis=1)
    if verbose:
        print('q=%d L=%d |A|=%d point orbits=%d monomial orbits=%d (%.1fs)' % (q, L, N, len(preps), len(mreps), time.time() - t0), flush=True)
    w = root_of_unity(q, p)
    # value tables for each shift g: V_g[point_rep, mon_rep] in {-1,0,1}
    Vg = []
    for g in range(q):
        sh = [shift(k, g, q) for k, _ in preps]
        nz, ng = masks(sh, q)
        zero = (asupp[None, :] & ~nz[:, None]) != 0
        par = popcount_parity((a1[None, :] & ng[:, None]))
        val = np.where(zero, 0, 1 - 2 * par).astype(np.int64)
        Vg.append(val)
    Htot = {}
    for s in range(q):
        rows = [i for i, (k, sz) in enumerate(preps) if sz == q or s == 0]
        coef = [pow(w, (-s * g) % q, p) for g in range(q)]
        M = np.zeros((len(rows), len(mreps)), dtype=object)
        acc = np.zeros((len(rows), len(mreps)), dtype=np.int64)
        for g in range(q):
            acc = (acc + Vg[g][rows, :] * (coef[g] % p)) % p   # coef < p < 2^59, |val|<=1 -> fits
        flat = [int(x) for x in acc.ravel()]
        Mm = flint.nmod_mat(len(rows), len(mreps), flat, p)
        R, rank = Mm.rref()
        # pivot columns
        Hs = {}
        r = 0
        ncol = len(mreps)
        # find pivots row by row
        Rn = np.array([[int(R[i, j]) for j in range(ncol)] for i in range(rank)], dtype=object) if rank * ncol < 4e7 else None
        piv = []
        for i in range(rank):
            row = Rn[i] if Rn is not None else [int(R[i, j]) for j in range(ncol)]
            for j in range(ncol):
                if row[j] != 0:
                    piv.append(j)
                    break
        for j in piv:
            Hs[int(mdeg[j])] = Hs.get(int(mdeg[j]), 0) + 1
        for dg, c in Hs.items():
            Htot[dg] = Htot.get(dg, 0) + c
        if verbose:
            print('  s=%d rows=%d rank=%d maxdeg=%d (%.1fs)' % (s, len(rows), rank, max(Hs) if Hs else -1, time.time() - t0), flush=True)
    H = [Htot.get(i, 0) for i in range(max(Htot) + 1)]
    return N, H


if __name__ == '__main__':
    q = int(sys.argv[1]); L = int(sys.argv[2])
    QR = sorted({(x * x) % q for x in range(1, q)})
    D = [0] + QR
    p = find_prime(q, 1 << 58)
    N, H = hilbert_cyclic(q, D, L, p)
    print('RESULT q=%d D=%s L=%d p=%d |A|=%d sumH=%d rho=%d H=%s' % (q, D, L, p, N, sum(H), len(H) - 1, H))
