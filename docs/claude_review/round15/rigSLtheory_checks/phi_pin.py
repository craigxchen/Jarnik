"""Round 15 (rigSLtheory), Theorem B: the label-independent angle set Phi(q0).

One-modulus setting: genuine labels on F' = F \\ {q0} (all levels), anything at q0.  P0 = {v in V_K (v ~ -v) :
U_1(v) = U_3(v) = 1}, i.e. D_1(Z_v) D_3(Z_v) = Re(Z_v^2) is fully supported on F' (all levels).  By Theorem B every SL
angle phi satisfies sin 2<v,phi> = sin 2 x_0(v) for v in P0, whatever the label at q0 (a label colliding v at q0 with
an odd unit has an A0 collision and is infeasible).  With psi = 2 phi mod 2 pi:
    <v, psi> = 2 x_0(v)   or   pi - 2 x_0(v)     (mod 2 pi),   v in P0.
Pick r independent v_1..v_r in P0 (matrix B, d = |det B|).  Every solution is psi = B^{-1}(b(tau) + 2 pi k),
tau in {0,1}^r, k in Z^r / B Z^r: 2^r d candidates.  Each candidate is excluded RIGOROUSLY if some w in P0 has
sin<w,psi> - sin 2x_0(w) enclosure (mpmath.iv, 40 digits) not containing 0.  psi_0 = 2 phi_0 (tau = 0, k in B Z^r)
always survives.  Output: number of candidates and the survivors.
Usage: python3 phi_pin.py FFILE q0 K p1..pr
"""
import sys, itertools
from fractions import Fraction as Fr
import numpy as np
from mpmath import iv
import sympy
from sllib import two_sq, Zv, Ds, VK, vq
from classify import load_F
iv.dps = 40


def P0_set(F, q0, K, ps):
    pis = [two_sq(p) for p in ps]; out = []
    for v in VK(len(ps), K):
        if next(x for x in v if x) < 0: continue
        z = Zv(pis, v); ok = True
        for s in (1, 3):
            d = abs(Ds(z, s)); g = 1
            for q in F:
                if q != q0: g *= q**vq(d, q)
            if g != d: ok = False; break
        if ok: out.append(v)
    return out


def main():
    Fpath, q0, K = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]); ps = [int(x) for x in sys.argv[4:]]
    F = load_F(Fpath); r = len(ps)
    P0 = P0_set(F, q0, K, ps)
    if sympy.Matrix(P0).rank() < r if P0 else True:
        print(f"ps={ps} K={K} q0={q0}: |P0|={len(P0)}, rank < r: Phi is not finite"); return
    cand = sorted(P0, key=lambda v: sum(map(abs, v)))[:40]
    best = None
    for sub in itertools.combinations(range(len(cand)), r):
        d = abs(round(np.linalg.det(np.array([cand[i] for i in sub], dtype=float))))
        if d and (best is None or d < best[0]):
            best = (d, sub)
            if d == 1: break
    d, sub = best
    B = sympy.Matrix([cand[i] for i in sub]); assert abs(B.det()) == d
    Binv = B.inv()
    args = [iv.atan2(iv.mpf(b), iv.mpf(a)) for a, b in map(two_sq, ps)]
    def x0(v): return sum((int(c)*a for c, a in zip(v, args) if c), iv.mpf(0))
    # coset representatives of Z^r / B Z^r: distinct Binv k mod Z^r
    reps = {}
    for k in itertools.product(range(d), repeat=r):
        key = tuple((x - sympy.floor(x)) for x in (Binv*sympy.Matrix(k)))
        reps.setdefault(key, k)
    assert len(reps) == d, (len(reps), d)
    rows = [B.row(i) for i in range(r)]
    surv = []
    tests = sorted(P0, key=lambda v: sum(map(abs, v)))
    for tau in itertools.product((0, 1), repeat=r):
        b = [2*x0(list(rows[i])) + tau[i]*(iv.pi - 4*x0(list(rows[i]))) for i in range(r)]
        for k in reps.values():
            rhs = [b[i] + 2*iv.pi*int(k[i]) for i in range(r)]
            psi = [sum((iv.mpf(int(Binv[j, i].p))/int(Binv[j, i].q)*rhs[i] for i in range(r) if Binv[j, i] != 0), iv.mpf(0))
                   for j in range(r)]
            excluded = False
            for w in tests:
                val = iv.sin(sum((c*p for c, p in zip(w, psi) if c), iv.mpf(0))) - iv.sin(2*x0(w))
                if val.a > 0 or val.b < 0: excluded = True; break
            if not excluded: surv.append((tau, k))
    print(f"ps={ps} K={K} q0={q0}: |P0|={len(P0)}, basis {[tuple(cand[i]) for i in sub]} det {d}, "
          f"candidates {2**r*d}, survivors {len(surv)}: {surv[:6]}")


if __name__ == '__main__':
    main()
