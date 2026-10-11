"""Exact verifier for the normalised problem.

verify(M, T, om, kappa, rho0, rho1) decides, in exact arithmetic, whether
    sigma_lambda(om) >= rho0 + rho1 ||lambda||_1 / 2      for every relevant lambda
(lambda in L_T \\ {0}, lambda.k in Z, k = (0, kappa)), and returns the exact minimum of
sigma_lambda - r(lambda) over relevant lambda (with a minimiser), or the violators found.

Finite reduction (Lemma 3.2 of exact.md): with an exact left inverse A (A Q = I),
    ||lambda||_inf <= k_inf N(lambda),   ||lambda||_1 <= k_1 N(lambda),
so sigma_lambda < r(lambda) (or <=) forces N(lambda) <= N_max = (Omega/2 + rho0)/(1 - rho1 k_1/2)
(when rho1 k_1 < 2), hence ||lambda||_inf <= k_inf N_max and ||lambda||_1 <= k_1 N_max.  All lattice
points of every relevant class in that box are enumerated.  If rho1 k_1 >= 2 the certificate is
inconclusive and the function says so.
"""
import math
from fractions import Fraction
import numpy as np
from core import (Group, Qmat, frac_part, left_inverse_lp, exact_left_inverse, coercivity, enum_class,
                  lcm_list, rows_of)


def integer_kernel_vector(Q, n=None):
    """A nonzero integer mu with Q mu = 0 (Q rank deficient), exact."""
    m = Q.shape[0]
    n = Q.shape[1] if n is None else n
    if m == 0:
        return [1] + [0] * (n - 1)
    # rational nullspace via elimination
    A = [[Fraction(int(Q[i, j])) for j in range(n)] for i in range(m)]
    piv = []
    r = 0
    for c in range(n):
        p = next((i for i in range(r, m) if A[i][c] != 0), None)
        if p is None:
            continue
        A[r], A[p] = A[p], A[r]
        pv = A[r][c]
        A[r] = [a / pv for a in A[r]]
        for i in range(m):
            if i != r and A[i][c] != 0:
                f = A[i][c]
                A[i] = [a - f * b for a, b in zip(A[i], A[r])]
        piv.append(c)
        r += 1
    free = [c for c in range(n) if c not in piv]
    fc = free[0]
    mu = [Fraction(0)] * n
    mu[fc] = Fraction(1)
    for i, c in enumerate(piv):
        mu[c] = -A[i][fc]
    L = lcm_list([v.denominator for v in mu])
    return [int(v * L) for v in mu]


def verify(M, T, om, kappa, rho0=0, rho1=0, strict=False, max_viol=50, A_float=None, verbose=False,
           max_group=200000):
    T = list(T)
    om = [Fraction(w) for w in om]
    rho0 = Fraction(rho0); rho1 = Fraction(rho1)
    kappa = [int(v) for v in kappa]
    Om = sum(om)
    Q = Qmat(M, T)
    out = {'M': M, 'nT': len(T), 'Omega': Om}
    G = Group(M, T)
    if not G.full_rank:
        mu = integer_kernel_vector(Q, M - 1)
        out.update(ok=False, reason='affinely dependent rows', violator=mu, min_slack=-Om / 2 - rho0)
        return out
    out['group'] = G.invs
    if G.order > max_group:
        out.update(ok=None, reason='group too large')
        return out
    if A_float is None:
        A_float = left_inverse_lp(Q, [float(w) for w in om])
    A = exact_left_inverse(Q, A_float)
    k_inf, k_1 = coercivity(Q, om, A)
    if rho1 > 0 and (rho1 * k_1 >= Fraction(3, 2) or k_1 > 1):
        from nu1 import kappa1_bound
        kb = kappa1_bound(M, T, om)
        if kb is not None and kb < k_1:
            k_1 = kb
    k_inf = min(k_inf, k_1 / 2)        # ||lambda||_inf <= ||lambda||_1 / 2 for zero-sum lambda
    out['k_inf'] = k_inf
    out['k_1'] = k_1
    if rho1 * k_1 >= 2:
        out.update(ok=None, reason='coercivity certificate inconclusive (rho1 k_1 >= 2)')
        return out
    Nmax = (Om / 2 + rho0) / (1 - rho1 * k_1 / 2)
    Binf = k_inf * Nmax
    B1 = k_1 * Nmax
    out['Binf'] = Binf
    out['B1'] = B1
    # integer scaling
    E = lcm_list([w.denominator for w in om])
    Wint = np.array([int(w * E) for w in om], dtype=np.int64)
    SW = int(Wint.sum())
    F0 = rho0.denominator; a0 = rho0.numerator
    F1 = rho1.denominator; a1 = rho1.numerator
    QT = Q.T.astype(np.int64)
    best = None
    viol = []
    nclasses = 0
    npoints = 0
    # classes: integral + kernel classes
    classes = [([0] * (M - 1), 1)]
    for nv in G.elements():
        if not any(nv):
            continue
        mu = G.element(nv)
        if sum(m_ * k_ for m_, k_ in zip(mu, kappa)) % 1 != 0:
            continue
        f = frac_part(mu)
        D = lcm_list([v.denominator for v in f])
        classes.append(([int(v * D) for v in f], D))
    out['relevant_rational_classes'] = len(classes) - 1
    for F, D in classes:
        nclasses += 1
        Binf_num = math.floor(Binf * D)
        B1_num = math.floor(B1 * D)
        scale = 2 * E * D * F0 * F1
        for V in enum_class(F, D, Binf_num, B1_num):
            npoints += V.shape[0]
            LS = np.abs(V @ QT)                          # D |lambda(S)|
            l1 = np.abs(V).sum(axis=1) + np.abs(V.sum(axis=1))   # D ||lambda||_1
            val = (2 * F0 * F1) * (LS @ Wint) - F0 * F1 * D * SW - 2 * E * D * a0 * F1 - E * a1 * F0 * l1
            # val = scale * (sigma - r)
            j = int(np.argmin(val))
            m = Fraction(int(val[j]), scale)
            if best is None or m < best[0]:
                lam = [Fraction(int(v), D) for v in V[j]]
                best = (m, lam)
            badmask = (val <= 0) if strict else (val < 0)
            if badmask.any() and len(viol) < max_viol:
                for jj in np.where(badmask)[0][: max_viol - len(viol)]:
                    viol.append((Fraction(int(val[jj]), scale), [Fraction(int(v), D) for v in V[jj]]))
    out['classes_checked'] = nclasses
    out['points_checked'] = npoints
    out['min_slack'] = best[0] if best else None
    out['argmin'] = best[1] if best else None
    out['violators'] = viol
    if best is None:
        out['ok'] = True
    else:
        out['ok'] = (best[0] > 0) if strict else (best[0] >= 0)
    return out


def lam_full(mu):
    """lambda (length M) from mu."""
    return [-sum(mu)] + list(mu)


if __name__ == '__main__':
    # smoke test: M=4 Hadamard core
    from core import pattern_of
    M = 4
    T = [pattern_of(M, r) for r in ([1, 2], [1, 3], [2, 3])]
    for kap in ([0, 0, 0], [1, 0, 0], [0, 0, 1]):
        r = verify(M, T, [1, 1, 1], kap, Fraction(1, 2), 0)
        print(kap, r['group'], r['relevant_rational_classes'], r['min_slack'], r['argmin'], r['ok'])
