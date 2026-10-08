"""Referee check 1 (Theorem C, step iii): the item-322/333 identities survive UNFLIPPED rows.

Model WF(t,b,F): rows x in F_2^t, b physical columns per nonzero label (integer weights w_j),
arbitrary orientations sigma_j, an arbitrary subset F of rows flipped (one flip per flipped row,
injective physical column). We check EXACTLY (Fractions / integers):

 (a) pair average:  E_x G_{x,x+d} = qhat(d) = sum_a (W_a - 4F_a/M) chi_a(d)   for all d != 0
 (b) item 322 (1):  V_c = W0/2 + (1/2) sum_{d in V\0} (-1)^{ell(d)} qhat(d)
                         + 2h Fcell/M - 2F/M + Fp - 2Fmatch
     with V_c = sum_j w_j |T_j|/2 computed from literal physical columns, for random
     subspaces V, functionals ell != 0 and cosets P -- including cosets containing unflipped rows
 (c) the composite of every coset character is a NONUNIT (some T_j != 0), so the phase inequality
     (4) is applicable on every coset (needed for the averaging (6)-(8) of item 333)
 (d) h*sum_{a|V=ell} q_a = Q - kappa + Psi  (kappa arbitrary), the identity behind (2a)
"""
import random, itertools
from fractions import Fraction as Fr

random.seed(2026)

def dot(a, x):
    return bin(a & x).count("1") & 1

def chi(a, x):
    return -1 if dot(a, x) else 1

def rand_subspace(t, s):
    while True:
        basis = [random.randrange(1, 2 ** t) for _ in range(s)]
        span = {0}
        for v in basis:
            span |= {u ^ v for u in span}
        if len(span) == 2 ** s:
            return basis, sorted(span)

def run(t, b, trials):
    M = 2 ** t
    labels = list(range(1, M))
    cols = []  # (label, sigma, weight)
    for a in labels:
        for _ in range(b):
            cols.append((a, random.choice([1, -1]), random.randint(1, 40)))
    ncol = len(cols)
    W0 = sum(c[2] for c in cols)
    W = {a: 0 for a in labels}
    for a, s, w in cols:
        W[a] += w
    cnt = 0
    for _ in range(trials):
        frac = random.choice([0.0, 0.1, 0.5, 0.9, 1.0])
        Fset = [x for x in range(M) if random.random() < frac]
        jcols = random.sample(range(ncol), len(Fset))  # injective flip columns
        flipcol = dict(zip(Fset, jcols))
        def s(x, j):
            a, sg, _ = cols[j]
            v = sg * chi(a, x)
            return -v if flipcol.get(x) == j else v
        f = {x: (cols[flipcol[x]][2] if x in flipcol else 0) for x in range(M)}
        alab = {x: cols[flipcol[x]][0] for x in flipcol}
        Fa = {a: 0 for a in labels}
        for x in flipcol:
            Fa[alab[x]] += f[x]
        Ftot = sum(Fa.values())
        q = {a: Fr(W[a]) - Fr(4 * Fa[a], M) for a in labels}
        qhat = lambda d: sum(q[a] * chi(a, d) for a in labels)
        S = [[s(x, j) for j in range(ncol)] for x in range(M)]
        G = lambda x, y: sum(cols[j][2] * S[x][j] * S[y][j] for j in range(ncol))
        # (a)
        for d in range(1, M):
            avg = Fr(sum(G(x, x ^ d) for x in range(M)), M)
            assert avg == qhat(d), ("pair average", t, d)
        # (b),(c),(d)
        for _r in range(6):
            sdim = random.randint(1, t)
            basis, V = rand_subspace(t, sdim)
            h = len(V)
            # random nonzero functional ell on V via values on basis
            while True:
                vals = [random.randint(0, 1) for _ in basis]
                if any(vals):
                    break
            ell = {}
            for coeffs in itertools.product([0, 1], repeat=sdim):
                v = 0; e = 0
                for c, bv, val in zip(coeffs, basis, vals):
                    if c:
                        v ^= bv; e ^= val
                ell[v] = e
            # labels restricting to ell on V: a.v = ell(v) for all v in V
            cell = [a for a in labels if all(dot(a, v) == ell[v] for v in V)]
            Fcell = sum(Fa[a] for a in cell)
            x0 = random.randrange(M)
            P = [x0 ^ v for v in V]
            c = {x0 ^ v: (-1 if ell[v] else 1) for v in V}
            T = [sum(c[x] * S[x][j] for x in P) for j in range(ncol)]
            Vc = Fr(sum(cols[j][2] * abs(T[j]) for j in range(ncol)), 2)
            Fp = sum(f[x] for x in P)
            Fmatch = sum(f[x] for x in P if x in alab and all(dot(alab[x], v) == ell[v] for v in V))
            rhs = Fr(W0, 2) + Fr(1, 2) * sum((-1) ** ell[d] * qhat(d) for d in V if d) \
                + Fr(2 * h * Fcell, M) - Fr(2 * Ftot, M) + Fp - 2 * Fmatch
            assert Vc == rhs, ("item-322 identity", t, frac)
            assert any(T), "unit composite"
            kappa = Fr(random.randint(-50, 50), 7)
            Psi = sum((-1) ** (ell[d] + 1) * (kappa - qhat(d)) for d in V if d)
            Q = sum(q.values())
            assert h * sum(q[a] for a in cell) == Q - kappa + Psi
            cnt += 1
    return cnt

tot = 0
for t, b, tr in [(3, 5, 25), (4, 5, 12), (5, 5, 4)]:
    tot += run(t, b, tr)
print("partial-flip identities (pair average, item-322 (1), nonunit composite, (2a) identity): %d coset checks PASS" % tot)
