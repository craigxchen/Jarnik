"""Round 15 (rigSLtheory): rigorous interval branch-and-bound for the slab system (1.2) of sharp.md at norm <= K,
one angle vector phi in [0, pi]^r, in the one-modulus setting (genuine labels on F \\ {q0} at all levels,
label ell at q0, level 1 there).

Constraint for (v, s), with y = 2<v, phi> (Lemma J0, exact rational box):
   2 m2^2/N - 1 <= cos y <= 1 - 2 m0^2/N,     m3^2/N - 1 <= sin y <= 1 - m1^2/N.
A box P of phi is discarded only if, for some v, the interval enclosure (mpmath.iv, outward rounded, 40 digits)
of cos y or sin y over P lies strictly outside the allowed interval: a rigorous proof that no phi in P works.
INFEASIBLE output = every box of [0, pi]^r discarded = certificate (interval arithmetic) that no phi exists.
(phi -> phi + pi e_j leaves every condition invariant, so [0, pi]^r suffices; v and -v give the same condition.)
If a box of width < tol has its centre satisfying every condition with rigorous interval evaluation, it is
reported as a verified witness phi.
Usage: python3 slab_iv.py FFILE q0 K "ell" p1 .. pr      (ell comma-separated level-1 label at q0)
"""
import sys, os, math
from fractions import Fraction as Fr
from mpmath import iv, mpf
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sllib import two_sq, genuine, nu, Zv, Ds, VK, vq
from classify import load_F

iv.dps = 40


def constraints(Fpath, q0, K, ell, ps):
    F = load_F(Fpath); r = len(ps); pis = [two_sq(p) for p in ps]
    m = q0-1 if q0 % 4 == 1 else q0+1; n = m//4
    out = []
    for v in VK(r, K):
        first = next(x for x in v if x)
        if first < 0: continue
        z = Zv(pis, v); N = z[0]**2 + z[1]**2
        ms = []
        for s in range(4):
            d = Ds(z, s); g = 1
            for q in F:
                if q != q0: g *= q**vq(d, q)
            ms.append(g)
        dd = sum(a*b for a, b in zip(v, ell)) % m
        if dd % n == 0: ms[dd//n] *= q0
        box = (Fr(2*ms[2]**2, N)-1, 1-Fr(2*ms[0]**2, N), Fr(ms[3]**2, N)-1, 1-Fr(ms[1]**2, N))
        out.append((v, box))
    return out


def ivq(x):
    return iv.mpf(x.numerator)/x.denominator


def check_box(cons, lo, hi):
    """lo, hi: lists of Fractions (multiples of 1, phi = pi * t). returns index of a violated constraint or None"""
    P = [iv.pi*iv.mpf([mpf(a.numerator)/a.denominator, mpf(b.numerator)/b.denominator]) for a, b in zip(lo, hi)]
    # the endpoints a/b are converted by mpf at 40 digits; enlarge by 1e-35 to stay outward
    P = [x + iv.mpf([-1e-35, 1e-35]) for x in P]
    for idx, (v, (lC, hC, lS, hS)) in enumerate(cons):
        y = 2*sum((c*x for c, x in zip(v, P) if c), iv.mpf(0))
        if y.delta > 7: continue
        C = iv.cos(y); S = iv.sin(y)
        if C.a > ivq(hC).b or C.b < ivq(lC).a or S.a > ivq(hS).b or S.b < ivq(lS).a:
            return idx
    return None


def bnb(cons, r, maxnodes=200000, tol=Fr(1, 2**20)):
    stack = [([Fr(0)]*r, [Fr(1)]*r)]
    nodes = 0
    while stack:
        lo, hi = stack.pop(); nodes += 1
        if nodes > maxnodes: return 'ABORT', nodes, None
        if check_box(cons, lo, hi) is not None: continue
        w = [b-a for a, b in zip(lo, hi)]
        j = max(range(r), key=lambda i: w[i])
        if w[j] < tol:
            mid = [(a+b)/2 for a, b in zip(lo, hi)]
            if check_point(cons, mid): return 'FEASIBLE', nodes, mid
            return 'UNDECIDED', nodes, mid  # a tiny box survived without a verified point: no certificate
        mid = (lo[j]+hi[j])/2
        a = (lo, hi[:j]+[mid]+hi[j+1:]); b = (lo[:j]+[mid]+lo[j+1:], hi)
        stack.append(a); stack.append(b)
        if nodes % 5000 == 0: print(f"    nodes {nodes} stack {len(stack)}", file=sys.stderr, flush=True)
    return 'INFEASIBLE', nodes, None


def check_point(cons, t):
    P = [iv.pi*(iv.mpf(x.numerator)/x.denominator) for x in t]
    for v, (lC, hC, lS, hS) in cons:
        y = 2*sum((c*x for c, x in zip(v, P) if c), iv.mpf(0))
        C = iv.cos(y); S = iv.sin(y)
        if not (C.b <= ivq(hC).a and C.a >= ivq(lC).b and S.b <= ivq(hS).a and S.a >= ivq(lS).b):
            return False
    return True


if __name__ == '__main__':
    Fpath, q0, K = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    ell = [int(x) for x in sys.argv[4].split(',')]; ps = [int(x) for x in sys.argv[5:]]
    cons = constraints(Fpath, q0, K, ell, ps)
    # sanity: genuine angles must be feasible for the genuine label (checked by caller with ell = genuine)
    res, nodes, wit = bnb(cons, len(ps))
    print(f"ps={ps} K={K} q0={q0} ell={ell} constraints={len(cons)}: {res} after {nodes} nodes; witness t (phi=pi t) = {wit}")
