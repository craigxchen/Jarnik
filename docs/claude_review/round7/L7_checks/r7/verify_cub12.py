"""EXACT verification of a cub12 member (class: cub13 without plane 456; D23=2, D13=1, D12=1)
derived from verify_cub13.py.  Original docstring follows.
EXACT verification of a cub13 member over a number field K = Q(alpha).
usage: python3 verify_cub13.py minpoly_<...>.json
The point: (r25, r34, r35) rational from the json, w = alpha (root of the json polynomial).
Checks: residual = 0; contact structure (all 56 boundary divisors) via gcds in K[s];
smoothness of W at the point (gradient of residual != 0); dominance of P^1 x W -> Mbar_{0,7}
(exact rank 4 of the differential at a rational parameter s0)."""
import sys, json, itertools
from fractions import Fraction as Fr
from nf import Field, KE, AD, pgcd, pdivmod, pderiv, pdeg, psub, pmul, ptrim, peval
import cub13exact as C

if sys.argv[1].endswith('.json'):
    J = json.load(open(sys.argv[1]))
    pt = [Fr(x) for x in J['point']]
    coeffs = [Fr(c) for c in reversed(J['coeffs_leading_first'])]
    F = Field(coeffs)
    alpha = F.gen()
    print('K = Q(alpha), alpha root of degree-%d polynomial' % F.n, J['coeffs_leading_first'])
else:   # four rationals r25 r34 r35 w : K = Q
    pt = [Fr(x) for x in sys.argv[1:4]]
    F = Field([Fr(0), Fr(1)])
    alpha = KE(F, Fr(sys.argv[4]))
    print('K = Q; point (r25, r34, r35, w) =', [str(x) for x in pt], sys.argv[4])
r25, r34, r35 = [KE(F, x) for x in pt]
one = KE(F, 1)
d = C.build(r25, r34, r35, alpha, one)
print('[1] residual (nonzero for a general cub12 member):', not d['residual'].iszero())

# ---- contact analysis -------------------------------------------------------------------------
MOB = (Fr(3), Fr(-1), Fr(2), Fr(5))     # t = (a s + b)/(g s + e): generic rational Moebius
def remob(p):
    a, b, g, e = MOB
    out = [KE(F, 0)] * 4
    for k, ck in enumerate(p + [KE(F, 0)] * (4 - len(p))):
        # ck * (a s + b)^k (g s + e)^(3-k)
        q = [KE(F, 1)]
        for _ in range(k): q = pmul(q, [KE(F, b), KE(F, a)])
        for _ in range(3 - k): q = pmul(q, [KE(F, e), KE(F, g)])
        q = q + [KE(F, 0)] * (4 - len(q))
        out = [out[i] + ck * q[i] for i in range(4)]
    return ptrim(out)
cs = [remob([KE(F, x) if not isinstance(x, KE) else x for x in p]) for p in C.coords(d)]
cs.append([])  # label 6 at 0
lab = [1, 2, 3, 4, 5, 6]
P = {}
for i, j in itertools.combinations(range(6), 2):
    P[(i, j)] = P[(j, i)] = psub(cs[i], cs[j])
# gcd structure
def gcd_list(ps):
    g = ps[0]
    for q in ps[1:]: g = pgcd(g, q) if g else q
    return g
def strip(e, q):
    while True:
        g = pgcd(e, q)
        if pdeg(g) <= 0: return e
        e = pdivmod(e, g)[0]
expected = {}
for T in ['456', '356', '346', '256', '246', '236', '345', '245', '135', '134', '125', '124']:
    expected[frozenset(int(x) for x in T)] = 1
expected[frozenset({1, 6})] = 3
for T in ['14', '15', '13', '12']:
    expected[frozenset(int(x) for x in T)] = 1
expected[frozenset({2, 3})] = 2
E = {}
ok = True; total = 0
for size in range(2, 6):
    for T in itertools.combinations(range(6), size):
        g = gcd_list([P[(i, j)] for i, j in itertools.combinations(T, 2)])
        if pdeg(g) <= 0:
            e = []
        else:
            e = g
            for k in range(6):
                if k not in T: e = strip(e, P[(T[0], k)])
        dg = max(pdeg(e), 0)
        key = frozenset(x + 1 for x in T)
        E[key] = e
        if dg != expected.get(key, 0):
            ok = False; print('   class mismatch at', sorted(key), dg, expected.get(key, 0))
        if dg > 0:
            total += dg
            sq = pdeg(pgcd(e, pderiv(e))) <= 0
            tr = all(pdeg(pgcd(e, pderiv(P[(i, j)]))) <= 0 for i, j in itertools.combinations(T, 2))
            if not (sq and tr): ok = False; print('   non-transversal / non-reduced at', sorted(key))
keys = [k for k in E if pdeg(E[k]) > 0]
for k1, k2 in itertools.combinations(keys, 2):
    if pdeg(pgcd(E[k1], E[k2])) > 0: ok = False; print('   shared contact point', sorted(k1), sorted(k2))
# point s = infinity: leading coefficients
lead = [p[3] if len(p) == 4 else KE(F, 0) for p in cs]
inf_ok = all(not (lead[i] - lead[j]).iszero() for i, j in itertools.combinations(range(6), 2))
print('[2] contact structure: class == cub12 (12 planes once, D16 thrice, D23 twice, D12/D13/D14/D15 once, others 0):', ok)
print('    total contact points', total, '; all reduced, transversal, pairwise distinct; no contact at s=inf:', inf_ok)
print('    parametrisation reduced (gcd of c1..c5 = 1):', pdeg(gcd_list([p for p in cs[:5]])) <= 0)

# ---- smoothness and dominance via exact automatic differentiation ----------------------------
nv = 5  # variables: s, r25, r34, r35, w
def var(x, k):
    g = [KE(F, 0)] * nv; g[k] = KE(F, 1); return AD(x if isinstance(x, KE) else KE(F, x), g)
oneA = AD.const(F, 1, nv)
D = C.build(var(r25, 1), var(r34, 2), var(r35, 3), var(alpha, 4), oneA)
grad = D['residual'].g[1:]
print('[3] (not needed for cub12) gradient of residual (d/dr25, d/dr34, d/dr35, d/dw) nonzero:', not all(x.iszero() for x in grad))
s0 = Fr(2, 7)
S0 = var(KE(F, s0), 0)
vals = [C.ev(p, S0) for p in C.coords(D)]
aff = [vals[i] / vals[4] for i in range(4)]       # affine chart c5 = 1 (c5 = X5)
Jm = [[a.g[k] for k in range(nv)] for a in aff]   # 4 x 5 over K
# tangent space of P^1 x W: d/ds plus kernel of grad (3-dim): basis
gi = next(k for k in range(4) if not grad[k].iszero())
basis = [[KE(F, 1)] + [KE(F, 0)] * 4]
for k in range(4):
    if k == gi: continue
    v = [KE(F, 0)] * 5; v[1 + k] = KE(F, 1); v[1 + gi] = -grad[k] / grad[gi]
    basis.append(v)
M = [[sum((Jm[r][c] * b[c] for c in range(5)), KE(F, 0)) for b in basis] for r in range(4)]
def det(M):
    M = [row[:] for row in M]; n = len(M); dt = KE(F, 1)
    for c in range(n):
        p = next((r for r in range(c, n) if not M[r][c].iszero()), None)
        if p is None: return KE(F, 0)
        if p != c: M[c], M[p] = M[p], M[c]; dt = -dt
        dt = dt * M[c][c]; iv = M[c][c].inv()
        for r in range(c + 1, n):
            f = M[r][c] * iv
            M[r] = [M[r][k] - f * M[c][k] for k in range(n)]
    return dt
dd = det(M)
print('[4] differential of P^1 x W -> P^4 at s0 = %s has rank 4 (det != 0) [here: with residual-gradient constraint, informational]:' % s0, not dd.iszero())
