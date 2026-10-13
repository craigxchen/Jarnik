"""Genuine pair-collision data for actual Gaussian primes on a given 0/1 profile (exact).

For a profile a in {0,1}^(M x r) (rows a_x = (1+S_x)/2), Gaussian primes pi_j = A_j + i B_j of norm p_j,
the genuine points are zhat_x = prod_j pi_j^(a_xj) conj(pi_j)^(1-a_xj)  (windings k = 0).  These are
actual lattice points of x^2 + y^2 = N, N = prod p_j.  For an odd prime power Q = q^a (q not dividing N):
    zhat_x / zhat_y = prod_j t_j^(a_xj - a_yj),   t_j = pi_j / conj(pi_j) in T_Q,
so with g a generator of T_Q and ell_j = dlog t_j, kappa_x = <a_x, ell> mod m_Q:
    zhat_x = i^s zhat_y (mod Q)   iff   kappa_x - kappa_y = s * iota (mod m_Q),   i = g^iota.
Unit-1 collision: kappa_x = kappa_y.  Any-unit collision: kappa_x = kappa_y mod m_Q/4.
At most one unit collides per Q (units are distinct mod an odd q), and collisions are nested in a.

The collision modulus of a pair (unit 1) is m_xy = prod over colliding Q <= X of q (one factor q per
colliding level), an exact integer; the pair demand is D_xy = 2 log m_xy.
"""
import math
import numpy as np
from gauss7 import TGroup, odd_prime_powers_upto, mQ

_TG_CACHE = {}


def tgroup(q, a):
    key = (q, a)
    if key not in _TG_CACHE:
        _TG_CACHE[key] = TGroup(q, a)
    return _TG_CACHE[key]


def point_classes(Aprof, gA, gB, ps, Xmax):
    """Return list of (q, a, Q, m, iota, kappa) for every odd prime power Q <= Xmax.
    Aprof: 0/1 int array (M x r).  gA, gB, ps: lists (Gaussian primes and norms), all p > Xmax."""
    assert min(ps) > Xmax
    out = []
    for (q, a, Q) in odd_prime_powers_upto(Xmax):
        G = tgroup(q, a)
        ell = G.dlog_t(gA, gB, ps)
        kappa = (Aprof @ ell) % G.m
        out.append((q, a, Q, G.m, G.iota, kappa))
    return out


def collision_moduli(classes, X, M):
    """Exact integer collision moduli for all pairs, from point classes.
    Returns (m1, many, unitmat): m1[x][y] = unit-1 modulus, many = any-unit modulus (python ints in
    object arrays), unitmat[x,y,Qindex] in {-1, 0,1,2,3}: the unit exponent s of the collision at Q
    (-1 if none)."""
    m1 = np.ones((M, M), dtype=object)
    many = np.ones((M, M), dtype=object)
    for (q, a, Q, m, iota, kappa) in classes:
        if Q > X:
            continue
        diff = (kappa[:, None] - kappa[None, :]) % m
        c1 = diff == 0
        ca = (diff % (m // 4)) == 0
        np.fill_diagonal(c1, False)
        np.fill_diagonal(ca, False)
        m1[c1] = m1[c1] * q
        many[ca] = many[ca] * q
    return m1, many


def demand_matrix(mod):
    M = mod.shape[0]
    D = np.zeros((M, M))
    for x in range(M):
        for y in range(M):
            if x != y:
                D[x, y] = 2.0 * math.log(int(mod[x, y]))
    return D


def window_counts(classes, M, lo, hi, anyunit=False):
    """c_xy = number of PRIMES q in (lo, hi] (level 1) with a collision of (x,y)."""
    C = np.zeros((M, M), dtype=np.int64)
    for (q, a, Q, m, iota, kappa) in classes:
        if a != 1 or not (lo < q <= hi):
            continue
        diff = (kappa[:, None] - kappa[None, :]) % m
        c = (diff % (m // 4) == 0) if anyunit else (diff == 0)
        np.fill_diagonal(c, False)
        C += c
    return C


def stats(D):
    M = D.shape[0]
    iu = np.triu_indices(M, 1)
    vals = D[iu]
    lam = float(np.linalg.eigvalsh(D)[-1])
    return dict(mean=float(vals.mean()), max=float(vals.max()), lam=lam,
                maxrow=float(D.sum(1).max()), meanrow=float(D.sum(1).mean()))


def expected_mean_independent(M, X, anyunit=False):
    """Exact expected pair demand under independent uniform classes (no parity):
    E D = 2 sum_{Q <= X} log q * u/m_Q with u = 1 (unit 1) or 4 (any unit)."""
    u = 4 if anyunit else 1
    return sum(2 * math.log(q) * u / mQ(q, a) for (q, a, Q) in odd_prime_powers_upto(X))
