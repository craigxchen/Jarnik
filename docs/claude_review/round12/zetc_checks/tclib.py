"""Core library for round-12 ZE-TC computations on T = 1 slab boxes.

A point set A subset {-1,0,1}^d (coordinates k' = (k_1..k_d)), cut by slabs |k(S)| <= L,
k(S) = sum_{x in S} k_x with S a subset of range(d) (a sign column with pattern S, base point
outside S; the base-point pattern is S = range(d)).

nu_A(u) = max{ i : <k,u>^i restricted to A is NOT a polynomial of degree < i on A }.
TC(A,m) = { u : nu_A(u) < m }.

Exact mode: fmpq / fmpz arithmetic (python-flint).  Modular mode: nmod_mat with a large prime
(used only as evidence or as a fast filter; every claim in the write-up that rests on a computation
is re-run exactly unless labelled otherwise).
"""
import itertools
from fractions import Fraction
import flint

P1 = 2305843009213693951  # 2^61 - 1 (prime)
P2 = 4611686018427387847  # prime < 2^62


def box_points(d, cuts):
    """All k in {-1,0,1}^d with |sum_{x in S} k_x| <= L for every (S, L) in cuts."""
    pts = []
    for k in itertools.product((-1, 0, 1), repeat=d):
        ok = True
        for S, L in cuts:
            if abs(sum(k[x] for x in S)) > L:
                ok = False
                break
        if ok:
            pts.append(k)
    return pts


def cube_monomials(d, maxdeg=None):
    """Exponent vectors in {0,1,2}^d sorted by total degree (k^3 = k on the cube)."""
    mons = list(itertools.product((0, 1, 2), repeat=d))
    if maxdeg is not None:
        mons = [a for a in mons if sum(a) <= maxdeg]
    mons.sort(key=lambda a: (sum(a), a))
    return mons


def mon_eval(k, a):
    v = 1
    for kx, ax in zip(k, a):
        if ax:
            v *= kx ** ax
    return v


def linval(k, u):
    return sum(kx * ux for kx, ux in zip(k, u))


def nu_mod(pts, u, p=P1, mons=None):
    """nu_A(u) modulo p (u integer vector).  Returns (nu, rho) where rho = top degree of S/J.
    Single rref with columns [deg0 mons, f_1, deg1 mons, f_2, deg2 mons, ...]; f_i = l_u^i on A."""
    d = len(pts[0])
    if mons is None:
        mons = cube_monomials(d)
    lv = [linval(k, u) % p for k in pts]
    maxdeg = sum(mons[-1])
    cols = []  # list of (tag, column)
    bydeg = {}
    for a in mons:
        bydeg.setdefault(sum(a), []).append(a)
    for i in range(0, maxdeg + 2):
        if i >= 1:
            cols.append(('f', i, [pow(x, i, p) for x in lv]))
        for a in bydeg.get(i, []):
            cols.append(('m', i, [mon_eval(k, a) % p for k in pts]))
    M = flint.nmod_mat(len(pts), len(cols), [c[2][r] for r in range(len(pts)) for c in cols], p)
    R, rank = M.rref()
    # pivot columns: first nonzero in each row
    piv = []
    ncols = len(cols)
    for r in range(rank):
        for j in range(ncols):
            if int(R[r, j]) != 0:
                piv.append(j)
                break
    nu = 0
    rho = 0
    for j in piv:
        tag, i, _ = cols[j]
        if tag == 'f':
            nu = max(nu, i)
        else:
            rho = max(rho, i)
    return nu, rho


def nu_exact(pts, u, mons=None):
    """Exact version over Q (fmpq_mat rref).  u may contain Fractions/ints."""
    d = len(pts[0])
    if mons is None:
        mons = cube_monomials(d)
    lv = [sum(Fraction(kx) * Fraction(ux) for kx, ux in zip(k, u)) for k in pts]
    maxdeg = sum(mons[-1])
    bydeg = {}
    for a in mons:
        bydeg.setdefault(sum(a), []).append(a)
    cols = []
    for i in range(0, maxdeg + 2):
        if i >= 1:
            cols.append(('f', i, [x ** i for x in lv]))
        for a in bydeg.get(i, []):
            cols.append(('m', i, [Fraction(mon_eval(k, a)) for k in pts]))
    # clear denominators per column to use fmpz_mat rank-profile via fmpq_mat rref
    M = flint.fmpq_mat(len(pts), len(cols), [flint.fmpq(c[2][r].numerator, c[2][r].denominator)
                                              for r in range(len(pts)) for c in cols])
    R, rank = M.rref()
    nu = 0
    rho = 0
    for r in range(rank):
        for j in range(len(cols)):
            if R[r, j] != 0:
                tag, i, _ = cols[j]
                if tag == 'f':
                    nu = max(nu, i)
                else:
                    rho = max(rho, i)
                break
    return nu, rho


def hilbert(pts, p=P1, mons=None):
    """dim (S/J)_i for all i (mod p rank profile; equals the Q-value unless p divides a minor)."""
    d = len(pts[0])
    if mons is None:
        mons = cube_monomials(d)
    M = flint.nmod_mat(len(pts), len(mons), [mon_eval(k, a) % p for k in pts for a in mons], p)
    R, rank = M.rref()
    h = {}
    for r in range(rank):
        for j in range(len(mons)):
            if int(R[r, j]) != 0:
                dg = sum(mons[j])
                h[dg] = h.get(dg, 0) + 1
                break
    return h
