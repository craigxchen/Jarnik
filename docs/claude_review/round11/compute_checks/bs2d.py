"""Complete exact/certified determination of Bs(A,m) subset (C^*)^2 for d = 2 configurations.
Output: curve components (exact gcd + factorisation over Q, binomial test), isolated points (resultant
projections, exact factorisation, certified acb root balls, certified exclusion of spurious pairs),
and the real-torus points near e.
"""
import itertools, math, random
import flint
from core import *

def laurent(A, c, ctx):
    mins = [min(k[i] for k in A) for i in range(2)]
    dct = {}
    for k, ck in zip(A, c):
        if ck:
            dct[(k[0] - mins[0], k[1] - mins[1])] = ck
    return ctx.from_dict(dct)

def gcd_all(polys):
    g = polys[0]
    for f in polys[1:]:
        g = g.gcd(f)
        if g.is_constant():
            break
    return g

def to_fmpz_poly_in(f, var):
    """f in fmpz_mpoly with one variable effectively -> fmpz_poly in that variable"""
    dct = f.to_dict()
    deg = max(e[var] for e in dct)
    co = [0] * (deg + 1)
    for e, c in dct.items():
        co[e[var]] += int(c)
    return flint.fmpz_poly(co)

def resultant_proj(Q, ctx, var, seed, tries=3):
    """gcd over several random pairs of Res_{other}(F,G) in Z[w_var]; contains the projection of V(Q)."""
    rnd = random.Random(seed)
    other = 1 - var
    R = None
    names = ctx.names()
    for _ in range(tries):
        F = sum((rnd.randint(-9, 9) * q for q in Q), ctx.from_dict({}))
        G = sum((rnd.randint(-9, 9) * q for q in Q), ctx.from_dict({}))
        r = F.resultant(G, names[other])
        rp = to_fmpz_poly_in(r, var) if not r.is_zero() else None
        if rp is None:
            continue
        R = rp if R is None else R.gcd(rp)
    return R

def analyse2d(A, m, label="", pr=print, near=0.25):
    basis = kernel_basis(A, m)
    check_kernel(A, m, basis)
    ctx = flint.fmpz_mpoly_ctx.get(('w', 2), 'lex')
    P = [laurent(A, c, ctx) for c in basis]
    G = gcd_all(P)
    res = dict(dimV=len(basis), sup=max(max(abs(x) for x in c) for c in basis))
    comps = []
    if not G.is_constant():
        for f, e in G.factor()[1]:
            dct = f.to_dict()
            binom = len(dct) == 2
            through_e = sum(int(c) for c in dct.values()) == 0
            comps.append((str(f), e, binom, through_e))
    res['curves'] = comps
    Q = [p // G for p in P] if not G.is_constant() else P
    # projections of the isolated points
    R1 = resultant_proj(Q, ctx, 0, 1)
    R2 = resultant_proj(Q, ctx, 1, 2)
    res['R1fac'] = [(str(f) if f.degree() <= 3 else "deg%d" % f.degree(), e) for f, e in R1.factor()[1]]
    res['R2fac'] = [(str(f) if f.degree() <= 3 else "deg%d" % f.degree(), e) for f, e in R2.factor()[1]]
    flint.ctx.prec = 200
    roots1 = [r for f, e in R1.factor()[1] for r, mu in f.complex_roots()]
    roots2 = [r for f, e in R2.factor()[1] for r, mu in f.complex_roots()]
    # candidate pairs, certified exclusion: some Q_c(a,b) ball excludes 0
    Qd = [q.to_dict() for q in Q]
    pts = []
    for a in roots1:
        if a == 0 or a.contains(0):
            continue
        for b in roots2:
            if b.contains(0):
                continue
            excluded = False
            for qd in Qd[:40]:
                s = flint.acb(0)
                for (i, j), c in qd.items():
                    s += int(c) * a ** i * b ** j
                if not s.contains(0):
                    excluded = True
                    break
            if not excluded:
                pts.append((a, b))
    res['isolated'] = pts
    # real-torus points among isolated candidates
    rt = []
    for a, b in pts:
        on = abs(a).contains(1) and abs(b).contains(1)
        if on:
            ua = float(a.arg().mid()) if hasattr(a, 'arg') else None
            ub = float(b.arg().mid())
            rt.append((ua, ub))
    res['isolated_real_torus'] = rt
    return res, basis, P, G, Q
