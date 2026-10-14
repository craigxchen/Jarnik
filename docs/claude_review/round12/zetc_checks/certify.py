"""Exact certification tools.

* hilbert_exact(pts): Hilbert function of A over Q (fmpz_mat rank profile in degree order; incremental by degree).
* zeta_certificate(pts, funcs): exact LP dual certificate  zeta(A) <= sum_f y_f width_A(f)  where y >= 0 rational
  with sum_f y_f |<v,f>| >= 1 for every primitive difference direction v of A (checked exactly with Fractions).
  Justification: for a lattice zonotope Z = k0 + sum [0,l_i] v_i in A, width_Z(f) = sum l_i |<v_i,f>| <= width_A(f)
  for every linear f; multiply by y_f and sum; each v_i is a positive integer multiple of a primitive direction of
  A - A, so sum_f y_f |<v_i,f>| >= 1, giving sum l_i <= sum_f y_f width_A(f).
"""
from fractions import Fraction
import numpy as np
import flint
from scipy.optimize import linprog
from tclib import cube_monomials, mon_eval
from zono import prim_dirs


def hilbert_exact(pts, mons=None):
    d = len(pts[0])
    if mons is None:
        mons = cube_monomials(d)
    bydeg = {}
    for a in mons:
        bydeg.setdefault(sum(a), []).append(a)
    H = {}
    prev = 0
    cols = []
    for dg in sorted(bydeg):
        cols += bydeg[dg]
        M = flint.fmpz_mat(len(pts), len(cols), [mon_eval(k, a) for k in pts for a in cols])
        r = M.rank()
        H[dg] = r - prev
        prev = r
        if r == len(pts):
            break
    return H


def widths(pts, funcs):
    P = np.array(pts)
    return [int((P @ np.array(f)).max() - (P @ np.array(f)).min()) for f in funcs]


def zeta_certificate(pts, funcs, denom=10**6):
    d = len(pts[0])
    dirs = prim_dirs(pts)
    fw = widths(pts, funcs)
    F = np.array(funcs, dtype=float)
    V = np.array(dirs, dtype=float)
    M = np.abs(V @ F.T)  # |V| x |F|
    # dual: min sum y_f w_f  s.t.  M y >= 1, y >= 0
    res = linprog(np.array(fw, dtype=float), A_ub=-M, b_ub=-np.ones(len(dirs)), bounds=(0, None), method='highs')
    y = res.x
    # rationalise and repair: scale up slightly then verify exactly
    best = None
    for scale in (Fraction(1), Fraction(1001, 1000), Fraction(101, 100)):
        yq = [Fraction(int(round(float(yi) * denom)), denom) * scale for yi in y]
        ok = True
        for v in dirs:
            s = sum(yq[j] * abs(sum(vi * fj for vi, fj in zip(v, funcs[j]))) for j in range(len(funcs)) if yq[j] != 0)
            if s < 1:
                ok = False
                break
        if ok:
            val = sum(yq[j] * fw[j] for j in range(len(funcs)))
            best = (val, yq)
            break
    return res.fun, best, len(dirs)
