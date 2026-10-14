"""Search for T = 1 slab-box models that are LP-supercritical, zonotope-free and have V(A,m) != 0.

Model: d coordinates, columns = patterns S_p (subsets of range(d)), slab half-widths L_p.
For each model:
  * A (exact enumeration), widths w_p = width_A(k(S_p)) (exact), h_p = 2 w_p;
  * Hilbert function H_A (mod-p rank profile; rho = top degree);
  * zeta upper bound by LP: max sum_v x_v s.t. sum_v x_v |<v,f>| <= width_A(f) for f in F
    (F = coordinates, column functionals, and all pair/triple sums of coordinates), v over primitive
    directions in A - A; solved by scipy (floating; candidates are re-certified exactly elsewhere);
  * a gadget lower bound for zeta is not needed for the faithful test (upper bound < m suffices).
Admissible m: max(2*max_p h_p, zeta_ub_floor + 1) <= m <= rho.
"""
import itertools, sys, math
import numpy as np
from scipy.optimize import linprog
from tclib import box_points, hilbert, cube_monomials


def widths(pts, funcs):
    P = np.array(pts)
    out = []
    for f in funcs:
        v = P @ np.array(f)
        out.append(int(v.max() - v.min()))
    return out


def primitive_dirs(pts):
    P = np.array(pts, dtype=np.int64)
    S = set()
    # differences of points, made primitive and sign-normalised
    for i in range(len(P)):
        D = P[i] - P[i + 1:]
        for v in D:
            g = 0
            for x in v:
                g = math.gcd(g, int(abs(x)))
            if g == 0:
                continue
            v = tuple(int(x) // g for x in v)
            # normalise sign: first nonzero positive
            for x in v:
                if x != 0:
                    if x < 0:
                        v = tuple(-y for y in v)
                    break
            S.add(v)
    return sorted(S)


def zeta_lp(pts, funcs, fw, dirs):
    if not dirs:
        return 0.0
    F = np.array(funcs, dtype=float)
    V = np.array(dirs, dtype=float)
    M = np.abs(V @ F.T).T  # |F| x |V|
    res = linprog(-np.ones(len(dirs)), A_ub=M, b_ub=np.array(fw, dtype=float), bounds=(0, None), method='highs')
    return -res.fun


def analyze(d, cols, Ls, extra_funcs=True, verbose=True, mons=None):
    cuts = [(S, L) for S, L in zip(cols, Ls)]
    pts = box_points(d, cuts)
    colf = [tuple(1 if x in S else 0 for x in range(d)) for S in cols]
    cw = widths(pts, colf)
    h = [2 * x for x in cw]
    H = hilbert(pts, mons=mons)
    rho = max(H)
    funcs = [tuple(1 if x == y else 0 for x in range(d)) for y in range(d)] + colf
    if extra_funcs:
        for x, y in itertools.combinations(range(d), 2):
            f = [0] * d; f[x] = 1; f[y] = 1; funcs.append(tuple(f))
            f = [0] * d; f[x] = 1; f[y] = -1; funcs.append(tuple(f))
    fw = widths(pts, funcs)
    dirs = primitive_dirs(pts)
    z = zeta_lp(pts, funcs, fw, dirs)
    mlo = max(2 * max(h), int(math.floor(z + 1e-9)) + 1)
    res = dict(d=d, N=len(pts), widths=cw, rho=rho, H=[H.get(i, 0) for i in range(rho + 1)], zeta_ub=z, mlo=mlo)
    res['ok'] = mlo <= rho
    if res['ok']:
        res['dimV'] = {m: sum(H.get(i, 0) for i in range(m, rho + 1)) for m in range(mlo, rho + 1)}
    return res


if __name__ == '__main__':
    pass
