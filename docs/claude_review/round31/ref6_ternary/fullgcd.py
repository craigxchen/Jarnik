"""Exact degree of the FULL gcd of specialised raw marks F_j(t,t^L), q=27.
Shared roots of t^L - v and t^L - v' are roots of v - v' in V\\{0}, hence lie in F_{q^2} = F_729 (k <= 3).
So gcd = G_formal(t,t^L) * prod_{alpha in F_729} (t-alpha)^{extra_alpha}, with
extra_alpha = min_j sum_v e_j(v) ord_alpha(t^L - v) - sum_v min_j e_j(v) ord_alpha(t^L - v)."""
import numpy as np, itertools, math
import descent_check as dc
K = dc.K
F27 = dc.F27; I = dc.I
def hasse_ord(alpha, vco, L):
    """order of vanishing at alpha of f = t^L - sum_b vco[b] t^b (vco in F_729)"""
    # coefficients of f; Hasse derivative D^s f(alpha) = sum_j C(j,s) f_j alpha^(j-s)
    deg = max(L, len(vco)-1)
    f = [0]*(deg+1); f[L] = 1
    for b, c in enumerate(vco): f[b] = int(K.add[f[b], K.neg[c]])
    for s in range(deg+1):
        val = 0
        for j in range(s, deg+1):
            if f[j] == 0: continue
            cb = math.comb(j, s) % 3
            if cb == 0: continue
            term = K.mul(f[j], K.pow(alpha, j-s)) if alpha else (f[j] if j == s else 0)
            if cb == 2: term = K.neg[term]
            val = int(K.add[val, term])
        if val != 0: return s
    return deg
def extra(forms, signs, L):
    k = len(forms[0])
    vs = list(itertools.product(F27, repeat=k))
    E = np.array([[1 + e*dc.chi27[dc.lam_val(l, x)] for x in vs] for l, e in zip(forms, signs)])
    gmin = E.min(axis=0)
    vco = [[K.mul(I, xb) for xb in x] for x in vs]
    tot = 0; detail = {}
    for alpha in range(729):
        aL = K.pow(alpha, L) if alpha else 0
        # v(alpha)
        fib = []
        for idx, c in enumerate(vco):
            val = 0
            for b, cb in enumerate(c):
                if cb: val = int(K.add[val, K.mul(cb, K.pow(alpha, b)) if alpha else (cb if b == 0 else 0)])
            if val == aL: fib.append(idx)
        if len(fib) < 2: continue
        ords = {idx: hasse_ord(alpha, vco[idx], L) for idx in fib}
        full = min(sum(E[j, idx]*ords[idx] for idx in fib) for j in range(len(forms)))
        formal = sum(gmin[idx]*ords[idx] for idx in fib)
        if full > formal:
            detail[alpha] = full - formal; tot += full - formal
    return tot, detail, int(gmin.sum())
if __name__ == '__main__':
    import sys
    t, d, g = extra([(1,0),(0,1),(1,1),(1,2)], [1]*4, 11)
    print('k=2 four-mark L=11: formal g=%d, extra=%d at alphas %s' % (g, t, d))
