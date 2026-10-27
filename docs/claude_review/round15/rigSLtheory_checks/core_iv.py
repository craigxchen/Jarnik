"""Deletion filter: a minimal (w.r.t. single deletions) subsystem of the slab constraints that is still
rigorously INFEASIBLE (slab_iv.bnb).  Prints each kept v with its demands m_s, genuine D_s, cofactors U_s
(w.r.t. F minus q0) and the type of the q0-collision.  Usage: python3 core_iv.py FFILE q0 K ell p1..pr"""
import sys
from fractions import Fraction as Fr
from slab_iv import constraints, bnb
from classify import load_F
from sllib import two_sq, genuine, Zv, Ds, vq
Fpath, q0, K = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]); ell = [int(x) for x in sys.argv[4].split(',')]
ps = [int(x) for x in sys.argv[5:]]; r = len(ps)
cons = constraints(Fpath, q0, K, ell, ps)
assert bnb(cons, r)[0] == 'INFEASIBLE'
keep = list(range(len(cons)))
# delete in order of decreasing norm first
for i in sorted(range(len(cons)), key=lambda i: -sum(abs(x) for x in cons[i][0])):
    trial = [k for k in keep if k != i]
    if bnb([cons[k] for k in trial], r, maxnodes=20000)[0] == 'INFEASIBLE': keep = trial
print(f"ps={ps} K={K} q0={q0} ell={ell}: minimal infeasible core of {len(keep)} monomials (bnb nodes {bnb([cons[k] for k in keep], r)[1]})")
F = load_F(Fpath); pis = [two_sq(p) for p in ps]; m = q0-1 if q0 % 4 == 1 else q0+1; n = m//4
_, L = genuine(ps, q0)
for k in keep:
    v = cons[k][0]; z = Zv(pis, v); N = z[0]**2+z[1]**2
    D = [Ds(z, s) for s in range(4)]
    g = []
    for s in range(4):
        gg = 1
        for q in F:
            if q != q0: gg *= q**vq(D[s], q)
        g.append(gg)
    U = [abs(D[s])//g[s] for s in range(4)]
    dd = sum(a*b for a, b in zip(v, ell)) % m; gd = sum(a*b for a, b in zip(v, L)) % m
    su = dd//n if dd % n == 0 else None; gu = gd//n if gd % n == 0 else None
    typ = 'none' if su is None else ('GEN' if su == gu else ('C' if q0*g[su] <= abs(D[su]) else 'A'))
    print(f"  v={v} N={N} D={D} U={U} q0-unit(ell)={su} q0-unit(gen)={gu} type={typ}")
