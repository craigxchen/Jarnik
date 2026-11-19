"""Four-mark example: formal G(t,t^11) (direct product over F_729 of (t^11 - v)^{min_j e_j(v)}), its descent to F_9,
the extra part of the full gcd, and where it lives."""
import numpy as np, itertools, pickle
import descent_check as dc
from f9lib import F9P, divmod9, gcd9
K = dc.K; L = 11
forms = [(1,0),(0,1),(1,1),(1,2)]
poly = np.array([1], np.int64); gsum = 0
for x in itertools.product(dc.F27, repeat=2):
    e = min(1 + dc.chi27[dc.lam_val(l, x)] for l in forms)
    gsum += e
    vcoef = [K.mul(dc.I, xb) for xb in x]
    for _ in range(e):
        new = np.zeros(len(poly)+L, np.int64); new[L:] = poly
        for b, vb in enumerate(vcoef):
            if vb == 0: continue
            term = K.neg[K.mul(poly, vb)]
            new[b:b+len(poly)] = K.add[new[b:b+len(poly)], term]
        poly = new
G0 = dc.to_f9(poly); assert G0 is not None
print('formal g =', gsum, ' G(t,t^11) has F_9 coefficients, degree', G0.deg())
d = pickle.load(open('four_mark_L11.pkl', 'rb'))
G = F9P(*d['G'])
E, r = divmod9(G, G0); assert r.iszero()
print('extra factor degree', E.deg())
from flint import fmpz_mod_poly_ctx, fq_default_ctx, fq_default_poly_ctx
fp3 = fmpz_mod_poly_ctx(3); f9 = fq_default_ctx(modulus=fp3([1, 0, 1]), var='i'); R = fq_default_poly_ctx(f9)
ii = f9.gen()
Ef = R([int(a) + int(b)*ii for a, b in zip(E.re, E.im)])
c, fac = Ef.factor()
print('extra factorisation over F_9:', [(str(f), e) for f, e in fac])
