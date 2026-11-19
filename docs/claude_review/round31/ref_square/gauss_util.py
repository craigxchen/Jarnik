"""Own Gaussian-integer utilities for the referee checks (independent of Astra's scripts)."""
from math import gcd, isqrt
from itertools import combinations
from sympy import factorint, Matrix


def gmul(a, b):
    return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])


def gconj(a):
    return (a[0], -a[1])


def gnorm(a):
    return a[0]*a[0]+a[1]*a[1]


def gpow(a, k):
    r = (1, 0)
    for _ in range(k):
        r = gmul(r, a)
    return r


def gdivmod_exact(a, b):
    """a/b if exact in Z[i], else None."""
    n = gnorm(b)
    num = gmul(a, gconj(b))
    if num[0] % n or num[1] % n:
        return None
    return (num[0]//n, num[1]//n)


def gdiv_round(a, b):
    n = gnorm(b)
    num = gmul(a, gconj(b))
    q = ((2*num[0]+n)//(2*n), (2*num[1]+n)//(2*n))
    return q


def ggcd(a, b):
    while b != (0, 0):
        q = gdiv_round(a, b)
        qb = gmul(q, b)
        a, b = b, (a[0]-qb[0], a[1]-qb[1])
    return a


_pi_cache = {}


def split_prime(p):
    """Gaussian prime pi with N(pi)=p, p = 1 mod 4."""
    if p in _pi_cache:
        return _pi_cache[p]
    for x in range(1, isqrt(p)+1):
        y2 = p - x*x
        y = isqrt(y2)
        if y*y == y2 and y > 0:
            _pi_cache[p] = (x, y)
            return (x, y)
    raise ValueError(p)


UNITS = [(1, 0), (0, 1), (-1, 0), (0, -1)]


def profile(z, H):
    """For a primitive tuple z of norm H (all p | H split), return dict p -> (e_p, [a_ip]) and units.
    z_i = eps_i prod pi_p^{a_ip} conj(pi_p)^{e_p-a_ip}."""
    f = factorint(H)
    prof = {}
    for p, e in f.items():
        assert p % 4 == 1, (p, H)
        pi = split_prime(p)
        a = []
        for zi in z:
            k = 0
            w = zi
            while True:
                q = gdivmod_exact(w, pi)
                if q is None:
                    break
                w = q
                k += 1
            a.append(k)
        prof[p] = (e, a)
    units = []
    for i, zi in enumerate(z):
        prod = (1, 0)
        for p, (e, a) in prof.items():
            pi = split_prime(p)
            prod = gmul(prod, gmul(gpow(pi, a[i]), gpow(gconj(pi), e-a[i])))
        u = gdivmod_exact(zi, prod)
        assert u in UNITS, (zi, prod, u)
        units.append(UNITS.index(u))
    return prof, units


def columns(prof, anchor=4):
    cols = []
    for p, (e, a) in sorted(prof.items()):
        cols.append([a[j]-a[anchor] for j in range(5) if j != anchor])
    return cols


def lattice_index(cols):
    """gcd of 4x4 minors (0 if rank < 4), via Smith/hermite on integer matrix."""
    if len(cols) < 4:
        return 0
    M = Matrix(cols).T  # 4 x n
    if M.rank() < 4:
        return 0
    # gcd of 4x4 minors = product of invariant factors; use smith normal form
    from sympy.matrices.normalforms import smith_normal_form
    S = smith_normal_form(M)
    d = 1
    for i in range(4):
        d *= abs(S[i, i])
    return d


def rank_mod2(vectors, dim):
    rows = [sum(((v[k] & 1) << k) for k in range(dim)) for v in vectors]
    basis = []
    for r in rows:
        for b in basis:
            r = min(r, r ^ b)
        if r:
            basis.append(r)
    return len(basis)


def parity_dependent(prof, units, anchor=4, with_units=True):
    """Are the four anchored ratios dependent in U(Q)/U(Q)^2 ?
    Coordinates: exponents of g_p mod 2 and the torsion bit (eps ratio in {+-i})."""
    others = [j for j in range(5) if j != anchor]
    vecs = []
    ps = sorted(prof)
    for j in others:
        v = [(prof[p][1][j]-prof[p][1][anchor]) % 2 for p in ps]
        if with_units:
            v.append(((units[j]-units[anchor]) % 4) % 2)
        vecs.append(v)
    dim = len(ps) + (1 if with_units else 0)
    return rank_mod2(vecs, dim) < 4
