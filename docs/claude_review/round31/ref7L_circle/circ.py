"""Independent checker for (VII.4): 20h <= 9D + 18r, and the triple lemma 2h_T <= d_T + 2u_T.
Exact arithmetic over QQ(i) with sympy Poly; root locations numerically (mpmath, 60 digits) on the
squarefree parts of exact factorizations, so multiplicities are exact."""
from itertools import combinations
from functools import reduce
import sympy as sp, mpmath as mp
mp.mp.dps = 60
t = sp.Symbol('t')
K = sp.QQ_I
def P(e): return sp.Poly(e, t, domain=K)
def bar(f):
    return sp.Poly.from_list([sp.conjugate(c) for c in f.all_coeffs()], t, domain=K)
def contact(marks):
    D = marks[0].degree()
    return D - max((a-b).degree() for a, b in combinations(marks, 2))
def roots_mult(N):
    """list of (complex root, multiplicity) of N, exact multiplicities"""
    out = []
    _, fl = N.sqf_list()
    for f, m in fl:
        if f.degree() == 0: continue
        cs = [mp.mpc(mp.mpf(sp.Rational(sp.re(c)).p)/sp.Rational(sp.re(c)).q, mp.mpf(sp.Rational(sp.im(c)).p)/sp.Rational(sp.im(c)).q) for c in f.all_coeffs()]
        rs = mp.polyroots(cs, maxsteps=400, extraprec=400)
        for r in rs: out.append((r, m))
    return out
TOL = mp.mpf(10)**-25
def circles(rts):
    """candidate circles with real centre: return list of (centre, radius). Includes one through every root."""
    ups = [r for r, m in rts if mp.im(r) > 0]
    cand = []
    for a in ups:
        cand.append((mp.mpf(0), abs(a)))   # placeholder single-root circle (centre 0)
    for a, b in combinations(ups, 2):
        if abs(mp.re(a) - mp.re(b)) > TOL:
            c = (abs(a)**2 - abs(b)**2) / (2*(mp.re(a) - mp.re(b)))
            cand.append((c, abs(a - c)))
    return cand
def oncount(rts, c, rho):
    return sum(m for r, m in rts if abs(abs(r - c) - rho) < TOL*max(1, rho))
def best_r(rts):
    """min over circles with real centre of r = (#roots of N off circle)/2"""
    tot = sum(m for r, m in rts)
    best = max([oncount(rts, c, rho) for c, rho in circles(rts)] + [0])
    return (tot - best) // 2 if (tot - best) % 2 == 0 else (tot - best) / 2, best
def line_r(rts):
    """vertical lines Re t = const: best on-line count"""
    ups = [r for r, m in rts if mp.im(r) > 0]
    best = 0
    for a in ups:
        best = max(best, sum(m for r, m in rts if abs(mp.re(r) - mp.re(a)) < TOL))
    tot = sum(m for r, m in rts)
    return (tot - best) / 2
def reduce_gcd(marks):
    g = reduce(lambda a, b: a.gcd(b), marks)
    return [m.exquo(g) for m in marks], g
def analyze(marks, name='', triples=True, verbose=True):
    N = marks[0]*bar(marks[0])
    assert all(m*bar(m) == N for m in marks), 'norms differ'
    assert len({m.as_expr() for m in marks}) == len(marks), 'not distinct'
    marks, g = reduce_gcd(marks)
    assert g.degree() == 0
    N = marks[0]*bar(marks[0])
    D = marks[0].degree(); h = contact(marks)
    rts = roots_mult(N)
    assert sum(m for _, m in rts) == 2*D
    r, on = best_r(rts)
    rl = line_r(rts)
    res = dict(name=name, k=len(marks), D=D, h=h, r=r, rline=rl)
    if len(marks) == 5:
        res['slack'] = 9*D + 18*r - 20*h
    tri = []
    if triples:
        # triple lemma with the best circle for the five-tuple AND best circle per triple
        for T in combinations(range(len(marks)), 3):
            F, gT = reduce_gcd([marks[i] for i in T])
            d = F[0].degree(); hT = contact(F)
            NT = F[0]*bar(F[0]); rT = roots_mult(NT)
            u, _ = best_r(rT)
            tri.append((T, d, hT, u, d + 2*u - 2*hT))
        res['min_triple_slack'] = min(x[4] for x in tri)
        res['sum_d'] = sum(x[1] for x in tri)
    if verbose: print(res, flush=True)
    return res, tri
