"""1-D test (exact integer arithmetic). S finite subset of Z.
(a) gcd of {P_c : c in V(S,m)} via sympy;  (b) exhaustive search over products of cyclotomic polynomials
(= all integer polys with all roots roots of unity, up to sign/monomial) of degree <= span(S) supported in S."""
import sympy as sp, sys
from sympy import Poly
y = sp.symbols('y')
def Vbasis(S, m):
    M = sp.Matrix([[k**j for k in S] for j in range(m)])
    return M.nullspace()
def gcd_forms(S, m):
    B = Vbasis(S, m)
    polys = [Poly(sum(c[i]*y**(S[i]-min(S)) for i in range(len(S))), y) for c in B]
    g = polys[0]
    for p in polys[1:]:
        g = sp.gcd(g, p)
    return sp.factor(g.as_expr()), len(B)
def cyclo_products(maxdeg):
    cyc = []
    k = 1
    while k < 200:
        ph = Poly(sp.cyclotomic_poly(k, y), y)
        if ph.degree() <= maxdeg: cyc.append((k, ph))
        k += 1
    out = []
    def rec(i, P, deg):
        out.append(P)
        for t in range(i, len(cyc)):
            kk, ph = cyc[t]
            if deg + ph.degree() <= maxdeg:
                rec(t, P * ph, deg + ph.degree())
    rec(0, Poly(1, y), 0)
    return out
def zeta1d(S):
    S = sorted(S); s0 = S[0]; Sn = set(x - s0 for x in S); span = max(Sn)
    best = (0, None)
    for P in cyclo_products(span):
        supp = [mon[0] for mon in P.monoms()]
        lo = min(supp); hi = max(supp)
        for a in range(0, span - (hi - lo) + 1):
            if all((s - lo + a) in Sn for s in supp):
                # order at 1
                o = 0; Q = P
                while Q.eval(1) == 0:
                    Q = sp.quo(Q, Poly(y - 1, y)); o += 1
                if o > best[0]: best = (o, sp.factor(P.as_expr()))
                break
    return best
if __name__ == '__main__':
    S = [int(x) for x in sys.argv[1].split(',')]
    print('S =', S)
    for m in range(1, len(S)):
        g, dim = gcd_forms(S, m)
        print('  m', m, 'dim V(S,m)', dim, 'gcd of P_c:', g, flush=True)
    L = sorted(set(a - b for a in S for b in S) - {0})
    reach = max(len([k for k in S if k + l in S]) for l in L)
    print('  shift criterion certifies ZE^all(S,m) for m <=', min(reach, len(S) - 1))
    z = zeta1d(S)
    print('  zeta (max order of a cyclotomic-product form supported in S):', z[0], z[1], flush=True)
