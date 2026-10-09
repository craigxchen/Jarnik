"""The C_5 polynomial ansatz for balanced rational curves in Mbar_{0,6}.

Configuration: anchor point 0 fixed at y=0, moving points y_1..y_5 (cubic polynomials in s):
    y_i(s) = c_i (s - rho_{i-1,i}) (s - rho_{i,i+1}) (s - tau_i)        (indices mod 5)
Boundary realisation (each of the 25 boundary divisors met once):
    D_{0i}            : own root tau_i of y_i
    D_{0,i,i+1}       : common root rho_{i,i+1} of y_i, y_{i+1}   (5 cycle edges)
    D_{abc}, abc = [5] minus a diagonal {i,i+2}: triple collision y_a=y_b=y_c at s=sigma_abc
    D_{ij}            : the remaining (third) root of y_i - y_j
The triple conditions are linear in z_T (common values) given rho, sigma:
    row i:  sum_cyc (sigma_T2 - sigma_T3) z_T1 / A_{i,T1} = 0,  A_{iT} = (sigma_T - rho_{i-1,i})(sigma_T - rho_{i,i+1}).
"""
from fractions import Fraction as Fr
from polylib import *

LABELS = [1, 2, 3, 4, 5]
EDGES = [(1, 2), (2, 3), (3, 4), (4, 5), (1, 5)]
DIAG = [(1, 3), (2, 4), (3, 5), (4, 1), (5, 2)]
TRIPLES = [tuple(sorted(set(LABELS) - set(d))) for d in DIAG]   # (2,4,5),(1,3,5),(1,2,4),(2,3,5),(1,3,4)


def edge_key(i, j):
    return tuple(sorted((i, j)))


def nbr_edges(i):
    im = 5 if i == 1 else i - 1
    ip = 1 if i == 5 else i + 1
    return edge_key(im, i), edge_key(i, ip)


def triples_of(i):
    return [T for T in TRIPLES if i in T]


def matrix(rho, sigma):
    """rho: dict edge->value, sigma: dict triple->value. Returns 5x5 Fraction matrix rows i, cols T."""
    M = []
    for i in LABELS:
        e1, e2 = nbr_edges(i)
        Ts = triples_of(i)
        row = [Fr(0)] * 5
        for k in range(3):
            T1, T2, T3 = Ts[k], Ts[(k + 1) % 3], Ts[(k + 2) % 3]
            A = (sigma[T1] - rho[e1]) * (sigma[T1] - rho[e2])
            row[TRIPLES.index(T1)] = (sigma[T2] - sigma[T3]) / A
        M.append(row)
    return M


def solve_curve(rho, sigma):
    """Given (rho, sigma) with det = 0, return dict with c_i, tau_i, z_T and the cubics y_i."""
    M = matrix(rho, sigma)
    ns = nullspace_frac(M)
    if len(ns) != 1:
        return None
    z = dict(zip(TRIPLES, ns[0]))
    c = {}; tau = {}
    for i in LABELS:
        e1, e2 = nbr_edges(i)
        Ts = triples_of(i)
        A = {T: (sigma[T] - rho[e1]) * (sigma[T] - rho[e2]) for T in Ts}
        den = z[Ts[0]] / A[Ts[0]] - z[Ts[1]] / A[Ts[1]]
        if den == 0:
            return None
        kappa = (sigma[Ts[0]] - sigma[Ts[1]]) / den
        if kappa == 0:
            return None
        c[i] = 1 / kappa
        tau[i] = sigma[Ts[0]] - z[Ts[0]] * kappa / A[Ts[0]]
        # consistency with third triple
        t3 = sigma[Ts[2]] - z[Ts[2]] * kappa / A[Ts[2]]
        assert t3 == tau[i], 'inconsistent tau'
    y = {}
    for i in LABELS:
        e1, e2 = nbr_edges(i)
        y[i] = scal(c[i], prod([linpoly(1, -rho[e1]), linpoly(1, -rho[e2]), linpoly(1, -tau[i])]))
    return dict(c=c, tau=tau, z=z, y=y, rho=rho, sigma=sigma)
