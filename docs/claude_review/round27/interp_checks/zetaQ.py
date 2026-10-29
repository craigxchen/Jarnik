"""Exact, complete search for RATIONAL forms with no nondegenerate zero supported in a finite S subset Z^2:
P = +-Y^a prod_j Phi_{k_j}(Y^{v_j})  (v_j primitive, first nonzero coord > 0).  [Lemma: a Q-Laurent polynomial
all of whose zeros in G_m^2(C) are degenerate has this form, by Thm A + Galois closure.]
Completeness: the Newton polygon of the product is the Minkowski sum of the segments [0, deg(Phi_k) v]; its width in
each coordinate direction must be <= that of S, which bounds the search.  Order at e = number of factors with k=1."""
import sympy as sp, math, sys
y = sp.symbols('y')
def cyclo_coeffs(k):
    p = sp.Poly(sp.cyclotomic_poly(k, y), y)
    return [int(c) for c in reversed(p.all_coeffs())]  # coeff of y^i
def factor_poly(v, k):
    cs = cyclo_coeffs(k)
    return {(i * v[0], i * v[1]): c for i, c in enumerate(cs) if c != 0}
def mul(P, Q):
    R = {}
    for a, x in P.items():
        for b, y_ in Q.items():
            m = (a[0] + b[0], a[1] + b[1]); R[m] = R.get(m, 0) + x * y_
    return {m: c for m, c in R.items() if c != 0}
def zetaQ(S):
    Sset = set(S)
    xs = [s[0] for s in S]; ys = [s[1] for s in S]
    WX = max(xs) - min(xs); WY = max(ys) - min(ys)
    facs = []
    for k in range(1, 200):
        dk = sp.totient(k)
        if dk > max(WX, WY): continue
        for a in range(0, WX // dk + 1):
            for b in range(-(WY // dk), WY // dk + 1):
                if (a, b) == (0, 0) or (a == 0 and b < 0) or math.gcd(a, abs(b)) != 1: continue
                if a * dk > WX or abs(b) * dk > WY: continue
                facs.append(((a, b), k, factor_poly((a, b), k)))
    best = [0, None]
    def width(P):
        a = [m[0] for m in P]; b = [m[1] for m in P]
        return max(a) - min(a), max(b) - min(b)
    def fits(P):
        ms = list(P.keys()); m0 = ms[0]
        for s in Sset:
            a0 = (s[0] - m0[0], s[1] - m0[1])
            if all((m[0] + a0[0], m[1] + a0[1]) in Sset for m in ms): return True
        return False
    def dfs(P, start, order, fl):
        if order > best[0]: best[0] = order; best[1] = list(fl)
        for i in range(start, len(facs)):
            v, k, F = facs[i]
            Q = mul(P, F)
            wx, wy = width(Q)
            if wx > WX or wy > WY: continue
            if not fits(Q): continue
            dfs(Q, i, order + (1 if k == 1 else 0), fl + [(v, k)])
    dfs({(0, 0): 1}, 0, 0, [])
    return best
if __name__ == '__main__':
    import ast
    S = ast.literal_eval(sys.argv[1])
    print(zetaQ(S))
