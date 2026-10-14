"""Rigorous LP upper bound for zeta(A) (largest sum L_i of a lattice zonotope in A), T = 1 boxes.

Every zonotope k0 + sum_i [0,L_i] v_i in A satisfies, for every functional f with |f| <= b_f on A,
sum_i L_i |f(v_i)| <= 2 b_f.  If y_f >= 0 with sum_f y_f |f(v)| >= 1 for every nonzero v in A - A, then
sum_i L_i <= sum_f 2 b_f y_f.  Functionals: coordinates (b=1), base point sum (b=1), each slab (b=L).
We use one weight per class (coordinates / base / class of slabs); classes given by the caller.
Every v in {-2..2}^d with |v_x| <= 2, |sum v| <= 2, |v(P)| <= 2L is enumerated (superset of A - A); the distinct
feature vectors are collected, the LP solved, the dual rationalised and checked exactly on all feature vectors.
"""
import sys, itertools, time
from fractions import Fraction
import numpy as np
from scipy.optimize import linprog
from flib import *
from model import profile


def features(M, pats, L, classes):
    d = M - 1
    cols = [rel(P) for P in pats]
    Pm = np.zeros((len(cols), d), dtype=np.int64)
    for j, C in enumerate(cols):
        Pm[j, C] = 1
    feats = set()
    vals = np.arange(-2, 3, dtype=np.int64)
    # enumerate first 5 coords outer, rest vectorised
    d1 = min(d, 6)
    tail = np.array(list(itertools.product(range(-2, 3), repeat=d - d1)), dtype=np.int64) if d > d1 else np.zeros((1, 0), np.int64)
    for head in itertools.product(range(-2, 3), repeat=d1):
        V = np.concatenate([np.tile(np.array(head, np.int64), (len(tail), 1)), tail], axis=1)
        s = V.sum(axis=1)
        PV = V @ Pm.T
        ok = (np.abs(s) <= 2) & (np.abs(PV) <= 2 * L).all(axis=1) & (np.abs(V).sum(axis=1) > 0)
        if not ok.any():
            continue
        V = V[ok]; s = s[ok]; PV = PV[ok]
        f = [np.abs(V).sum(axis=1), np.abs(s)]
        for cl in classes:
            f.append(np.abs(PV[:, cl]).sum(axis=1))
        F = np.stack(f, axis=1)
        for row in np.unique(F, axis=0):
            feats.add(tuple(int(x) for x in row))
    return sorted(feats)


def bound(M, pats, L, classes):
    feats = features(M, pats, L, classes)
    caps = [2 * (M - 1), 2] + [2 * L * len(cl) for cl in classes]
    F = np.array(feats, dtype=float)
    res = linprog(np.array(caps, float), A_ub=-F, b_ub=-np.ones(len(F)), bounds=[(0, None)] * len(caps), method='highs')
    y = [Fraction(x).limit_denominator(10 ** 4) for x in res.x]
    worst = min(sum(Fraction(a) * b for a, b in zip(f, y)) for f in feats)
    assert worst > 0
    if worst < 1:
        y = [b / worst for b in y]
    assert all(sum(Fraction(a) * b for a, b in zip(f, y)) >= 1 for f in feats)
    return sum(Fraction(c) * b for c, b in zip(caps, y)), res.fun, y, len(feats)


if __name__ == '__main__':
    name = sys.argv[1]; L = int(sys.argv[2]); mult = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    M, pats = profile(name)
    q = M - 1
    # classes of slabs: orbits under the shift (patterns grouped by size as a proxy, checked to be shift-orbits)
    S = [tuple(sorted(x - 1 for x in P)) for P in pats]
    orbits = []
    seen = set()
    for i, P in enumerate(S):
        if i in seen:
            continue
        orb = set()
        cur = P
        for g in range(q):
            cur = tuple(sorted(((mult * x) % q if mult else (x + 1) % q) for x in cur))
            orb.add(S.index(cur))
        seen |= orb
        orbits.append(sorted(orb))
    t = time.time()
    b, f, y, nf = bound(M, pats, L, orbits)
    print('%s L=%d slab orbits=%s  LP bound zeta <= %s = %.4f (LP value %.4f), %d feature vectors, y=%s (%.1fs)' % (
        name, L, [len(o) for o in orbits], b, float(b), f, nf, y, time.time() - t))
