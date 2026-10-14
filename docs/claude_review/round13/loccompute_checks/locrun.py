"""Driver: J(A) versus local ideals for small models (d <= 8), mod-p sandwich (+ exact ranks on request)."""
import itertools, math, sys, time
import numpy as np
from loclib import Ring, hilbert_modp, hilbert_exact, image, jproj_gens, pullback, ideal_dims_modp, ideal_dims_exact, PR
import models


def width_dirs(A, d, R, maxdeg):
    """primitive v in [-R,R]^d (first nonzero > 0) with |image of <v,.> on A| <= maxdeg.  returns [(v, size)]"""
    A = np.array(A, dtype=np.int64)
    out = []
    vs = []
    for v in itertools.product(range(-R, R + 1), repeat=d):
        nz = [t for t in v if t]
        if not nz or nz[0] < 0:
            continue
        if math.gcd(*[abs(t) for t in nz]) != 1:
            continue
        vs.append(v)
    V = np.array(vs, dtype=np.int64)
    for s in range(0, len(V), 5000):
        blk = V[s:s + 5000]
        P = A @ blk.T
        P = np.sort(P, axis=0)
        sizes = 1 + (np.diff(P, axis=0) != 0).sum(axis=0)
        for t in np.nonzero(sizes <= maxdeg)[0]:
            out.append((tuple(int(z) for z in blk[t]), int(sizes[t])))
    return out


def width_gens(ring, wd, mod=PR):
    g = {}
    for v, s in wd:
        vec = ring.one(mod)
        for t in range(s):
            vec = ring.mul_lin(vec, t, v, mod)
        if mod is None:
            if any(int(z) != 0 for z in vec):
                g.setdefault(s, []).append(vec)
        elif vec.any():
            g.setdefault(s, []).append(vec)
    return g


def family_gens(ring, A, fams, maxdeg, mod=PR, cache=None, info=None):
    g = {}
    if cache is None:
        cache = {}
    for F in fams:
        P = image(A, F)
        if len(F) >= 3 and mod is not None:
            from jfast2 import jproj_gens_fast2
            jp, rho, _ = jproj_gens_fast2(P, len(F), maxdeg)
        else:
            jp, rho = jproj_gens(P, len(F), maxdeg)
        if info is not None:
            info.append((F, len(P), rho, {i: len(v) for i, v in jp.items() if v}))
        for i, forms in jp.items():
            for f in forms:
                vec = pullback(ring, F, f, i, mod, cache)
                if (mod is None and any(int(z) != 0 for z in vec)) or (mod is not None and vec.any()):
                    g.setdefault(i, []).append(vec)
    return g


def merge(*gs):
    out = {}
    for g in gs:
        for i, v in g.items():
            out.setdefault(i, []).extend(v)
    return out


def cum_eq(a, b):
    return [x == y for x, y in zip(itertools.accumulate(a), itertools.accumulate(b))]
