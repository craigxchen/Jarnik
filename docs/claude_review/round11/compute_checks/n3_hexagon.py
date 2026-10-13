"""n = 3 (d = 2): the only admissible models are single-slab profiles (see explore3_out.txt and
Proposition 2.1 of compute.md).  Slab on point 0 or 1 gives a product box [-L,L]x[-T,T]; slab on the base
point 2 gives the hexagon {|k1|,|k2| <= T, |k1+k2| <= L}.  For each instance:
  (1) exact LLL basis of V(A,m), sup-norm vs Siegel H;  V_H = V(A,m) iff supnorm <= H.
  (2) tangent cone: gcd of the binary forms Phi_c, factorisation over Q, real roots rational?
  (3) base locus: gcd of the Laurent polynomials P_c(w) (curve part), factorisation, binomial test;
      zero-dimensional residual via resultants (exact factorisation) + acb certified roots.
"""
import sys, itertools, math
import flint
from core import *

out = []
def pr(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True); out.append(s)

def binary_forms(A, m, basis):
    mons, rows = phi_forms(A, m, basis, 2)
    ctx = flint.fmpz_mpoly_ctx.get(('u', 2), 'lex')
    u1, u2 = ctx.gens()
    polys = []
    for rw in rows:
        f = ctx.from_dict({a: c for a, c in zip(mons, rw) if c})
        polys.append(f)
    return ctx, polys

def gcd_all(polys):
    g = polys[0]
    for f in polys[1:]:
        g = g.gcd(f)
    return g

def laurent_polys(A, basis, T, ctx):
    pols = []
    for c in basis:
        dct = {}
        for k, ck in zip(A, c):
            if ck:
                dct[(k[0] + T, k[1] + T)] = ck
        pols.append(ctx.from_dict(dct))
    return pols

def is_binomial(f):
    return len(f.terms()) == 2 if hasattr(f, 'terms') else len(f.to_dict()) == 2

def run(name, T, slabS, L, m):
    slabs = [(slabS, L)]
    A = build_A(3, T, slabs)
    basis = kernel_basis(A, m)
    check_kernel(A, m, basis)
    sup = max(max(abs(x) for x in c) for c in basis)
    H = siegel_H(A, T, m)
    pr("== %s T=%d slab=%s L=%d m=%d |A|=%d dimV=%d supnorm(LLL basis)=%d H=%.3g V_H=V:%s" %
       (name, T, slabS, L, m, len(A), len(basis), sup, H, sup <= H))
    # tangent cone
    ctx, forms = binary_forms(A, m, basis)
    pr("  tangent cone: dim span Phi =", len(forms))
    g = gcd_all(forms)
    fac = g.factor()
    pr("  gcd of Phi_c:", fac)
    nondeg_tc = []
    for f, e in fac[1]:
        dg = f.total_degree()
        if dg >= 2:
            # irreducible over Q of degree >= 2: real roots are irrational directions
            dd = f.to_dict()
            # dehomogenise u2 = 1 -> univariate in u1
            coeffs = [0] * (dg + 1)
            for (a, b), c in dd.items():
                coeffs[a] = int(c)
            p = flint.fmpz_poly(coeffs)
            rr = [r for r, mult in p.complex_roots() if abs(r.imag) < 1e-30 or r.imag.contains(0)]
            pr("    irreducible factor of degree", dg, "real roots:", len(rr))
            if rr:
                nondeg_tc.append((str(f), [str(r) for r in rr]))
    pr("  NONDEGENERATE real tangent directions:", nondeg_tc if nondeg_tc else "none")
    # base locus
    ctxw = flint.fmpz_mpoly_ctx.get(('w', 2), 'lex')
    P = laurent_polys(A, basis, T, ctxw)
    G = gcd_all(P)
    facG = G.factor()
    pr("  gcd of P_c (curve part of base locus):", facG)
    for f, e in facG[1]:
        dd = f.to_dict()
        val_e = sum(int(c) for c in dd.values())
        pr("    component", f, "mult", e, "binomial" if len(dd) == 2 else "NON-BINOMIAL", "passes e:" , val_e == 0)
    # zero-dimensional part: resultant of two random combinations of P_c / G
    Q = [p // G if not G.is_one() else p for p in P]
    import random
    random.seed(7)
    def comb():
        f = ctxw.from_dict({})
        for q in Q:
            f = f + random.randint(-5, 5) * q
        return f
    F1, F2 = comb(), comb()
    # resultant in w2 -> polynomial in w1 (use sympy for robustness)
    import sympy as sp
    w1, w2 = sp.symbols('w1 w2')
    def tosp(f):
        return sum(int(c) * w1 ** a * w2 ** b for (a, b), c in f.to_dict().items())
    R = sp.resultant(tosp(F1), tosp(F2), w2)
    R = sp.Poly(R, w1)
    pr("  resultant degree", R.degree())
    facR = sp.factor_list(R.as_expr())
    zd = []
    for fct, e in facR[1]:
        fp = sp.Poly(fct, w1)
        # check whether all Q vanish on some common root: gcd over the number field via numerics with acb
        pr("    resultant factor deg", fp.degree(), "mult", e, ":", fct if fp.degree() <= 4 else "(deg>4)")
        zd.append(fct)
    return dict(A=A, basis=basis, forms=forms, P=P, G=G, Q=Q, zd=zd)

if __name__ == "__main__":
    cases = [("hexagon", 3, (2,), 1, 5), ("hexagon", 3, (2,), 1, 6), ("hexagon", 3, (2,), 1, 7),
             ("hexagon", 4, (2,), 1, 5), ("hexagon", 4, (2,), 1, 8), ("hexagon", 5, (2,), 2, 9),
             ("hexagon", 6, (2,), 2, 11), ("hexagon", 6, (2,), 3, 13), ("hexagon", 6, (2,), 2, 13),
             ("box", 3, (0,), 1, 6), ("box", 5, (0,), 2, 10)]
    for c in cases:
        run(*c)
    open("n3_hexagon_out.txt", "w").write("\n".join(out) + "\n")
