"""Polynomial (function-field) two-class exponent.
Points y_i(t) in Q[t].  For each bipartition sigma, G_sigma = gcd of same-side differences (monic).
The exponent is E = (sum_sigma deg G_sigma + f_inf) / deg(diam), where deg(diam) = max_i,j deg(y_i - y_j)
and f_inf is the two-class order at infinity (homogeneous version; reported separately).
Lower bound transfer to integers: for coprime integer (alpha,beta) the homogenised values satisfy
q2(y(alpha,beta)) >= c * prod_sigma |G_sigma(alpha,beta)| (pairwise gcds of values of coprime forms
divide fixed resultants), so the integer ratio tends to E_aff = sum deg G_sigma / d along t -> infinity
(beta = 1), or to E_hom when both alpha, beta grow.
"""
import sympy as sp
from itertools import product
t = sp.symbols('t')

def split_polygcds(ys):
    B = len(ys); out = {}
    for m in range(1 << (B - 1)):
        sides = [0] + [(m >> (i - 1)) & 1 for i in range(1, B)]
        g = sp.Integer(0); rep = [None, None]
        for i in range(B):
            s = sides[i]
            if rep[s] is None: rep[s] = ys[i]
            else: g = sp.gcd(g, sp.expand(ys[i] - rep[s]))
        if g != 0:
            out[m] = sp.Poly(g, t)
    return out

def exponent(ys):
    ys = [sp.expand(y) for y in ys]
    d = max(sp.degree(sp.expand(a - b), t) for a in ys for b in ys if sp.expand(a - b) != 0)
    G = split_polygcds(ys)
    tot = sum(g.degree() for g in G.values())
    # order at infinity: max f such that <=2 classes "mod t^-f": same-side differences have degree <= d - f
    finf = 0
    B = len(ys)
    for m in range(1 << (B - 1)):
        sides = [0] + [(m >> (i - 1)) & 1 for i in range(1, B)]
        rep = [None, None]; dd = -1
        for i in range(B):
            s = sides[i]
            if rep[s] is None: rep[s] = ys[i]
            else:
                df = sp.expand(ys[i] - rep[s])
                dd = max(dd, sp.degree(df, t) if df != 0 else -1)
        if dd >= 0: finf = max(finf, d - dd)
    return d, tot, finf, {m: g.as_expr() for m, g in G.items() if g.degree() > 0}

if __name__ == "__main__":
    a_, b_ = t, sp.Integer(1)
    al, be = t, 1
    fam5 = [0, al * (al + be), be * (al + 2 * be), al * (al + be) + be * (al + 2 * be), (al + be) * (al + 2 * be)]
    print("B=5 family:", exponent(fam5))
