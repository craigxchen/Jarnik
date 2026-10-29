"""K1: 1-D example S={0,1,3,10,11,13} of interp Prop 3.3 -- exact recomputation (sympy), own code."""
import itertools, sys
import sympy as sp
from sympy import Matrix, symbols, Poly, gcd, factor_list, cyclotomic_poly, totient
y = symbols('y')
S = [0, 1, 3, 10, 11, 13]
out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)

def V(S, m):
    M = Matrix([[k ** t for k in S] for t in range(m)]) if m > 0 else Matrix.zeros(0, len(S))
    return M.nullspace() if m > 0 else [Matrix([1 if i == j else 0 for i in range(len(S))]) for j in range(len(S))]

for m in range(1, 6):
    B = V(S, m)
    polys = [sp.expand(sum(c * y ** k for c, k in zip(list(b), S))) for b in B]
    g = polys[0]
    for q in polys[1:]:
        g = gcd(g, q)
    g = sp.factor(g)
    pr("m=%d dim V=%d gcd=%s" % (m, len(B), g))
    if len(B) == 1:
        P = sp.Poly(sp.together(polys[0] * sp.lcm([sp.fraction(sp.nsimplify(c))[1] for c in list(B[0])])), y)
        cf = P.all_coeffs()
        pr("   P_c coefficients (high->low):", cf)
        fl = sp.factor_list(P.as_expr())
        pr("   factorisation:", fl)
        Q = sp.quo(P.as_expr(), (y - 1) ** 5)
        Qp = sp.Poly(Q, y)
        pr("   cofactor deg", Qp.degree(), "lead", Qp.LC(), "const", Qp.TC(), "|prod roots| =", abs(sp.Rational(Qp.TC(), Qp.LC())))
        rts = sp.Poly(Q, y).nroots(n=30)
        pr("   cofactor root moduli:", sorted([float(abs(r)) for r in rts]))
# (b) shift sizes
diffs = sorted(set(b - a for a in S for b in S if b != a))
cnt = {l: len([k for k in S if k + l in S]) for l in diffs}
pr("max |S_l| =", max(cnt.values()), "attained at", [l for l in cnt if cnt[l] == max(cnt.values())])
# (c) all rational products +- y^a prod Phi_k^{m_k} with support in S
ks = [k for k in range(1, 200) if totient(k) <= 13]
cyc = {k: sp.Poly(cyclotomic_poly(k, y), y) for k in ks}
phis = {k: int(totient(k)) for k in ks}
best = 0; bestex = []
def rec(i, deg, cur, mult1):
    global best, bestex
    # test current product
    co = cur.all_coeffs()[::-1]   # low -> high
    supp = [t for t, c in enumerate(co) if c != 0]
    for a in range(0, 14 - cur.degree()):
        if all((t + a) in S for t in supp):
            if mult1 > best:
                best = mult1; bestex = [(cur.as_expr(), a)]
            elif mult1 == best and len(bestex) < 5:
                bestex.append((cur.as_expr(), a))
            break
    for jj in range(i, len(ks)):
        k = ks[jj]
        if deg + phis[k] <= 13:
            rec(jj, deg + phis[k], cur * cyc[k], mult1 + (1 if k == 1 else 0))
rec(0, 0, sp.Poly(1, y), 0)
pr("max order at y=1 of a rational cyclotomic product supported in S:", best, "examples:", [(sp.factor(e), a) for e, a in bestex])
open("k1_onedim_out.txt", "w").write("\n".join(out) + "\n")
