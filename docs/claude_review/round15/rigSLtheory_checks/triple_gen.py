"""Theorem B(2) exact test: P0(q0) = {v in V_K : U_1(v) = U_3(v) = 1} (cofactors w.r.t. F \\ {q0}, all levels; P0 is
closed under v -> -v).  T = union of additive triples {v, w, v+w} in P0 (v+w != 0, v = w allowed).  Prints |P0|, |T|
and the index [Z^r : <T>] (gcd of r x r minors, via Smith normal form).  Index 1 => Phi(q0) = {phi_0} mod pi Z^r.
Usage: python3 triple_gen.py FFILE K QMAX p1..pr     (scans every q0 in F with q0 <= QMAX)"""
import sys
from sympy import Matrix
from sympy.matrices.normalforms import smith_normal_form
from sympy import ZZ
from sllib import two_sq, Zv, Ds, VK, vq
from classify import load_F
Fpath, K, QMAX = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]); ps = [int(x) for x in sys.argv[4:]]
F = load_F(Fpath); r = len(ps); pis = [two_sq(p) for p in ps]
V = VK(r, K)
DD = {v: (abs(Ds(Zv(pis, v), 1)), abs(Ds(Zv(pis, v), 3))) for v in V}
def index_of(vecs):
    if not vecs: return 0
    M = Matrix(vecs)
    if M.rank() < r: return 0
    S = smith_normal_form(M, domain=ZZ)
    prod = 1
    for i in range(r): prod *= abs(S[i, i])
    return prod
for q0 in [q for q in F if q <= QMAX]:
    Fp = [q for q in F if q != q0]
    P0 = set()
    for v in V:
        ok = True
        for d in DD[v]:
            g = 1
            for q in Fp:
                if d % q == 0: g *= q**vq(d, q)
            if g != d: ok = False; break
        if ok: P0.add(v)
    T = set()
    P0l = list(P0)
    for v in P0l:
        for w in P0l:
            u = tuple(a+b for a, b in zip(v, w))
            if any(u) and u in P0: T.update([v, w, u])
    idx = index_of(sorted(T))
    print(f"q0={q0}: |P0|={len(P0)} |T|={len(T)} index[Z^r:<T>]={idx if idx else 'rank<r'}", flush=True)
