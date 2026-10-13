# Leaf direction of F_G along pencils through m, as s -> 0 (R(0) = m). If the limit directions in P(T_m G)
# differ between pencils, the line field has no continuous extension at m, i.e. m ∈ Sing(F_G).
import sys, random; sys.path.insert(0, '.')
from reflib import *
import sympy as sp
cube = {(i,j,k): meet3(i,j,k) for i in (1,2) for j in (3,4) for k in (5,6)}
Ekeys = [(1,3,5),(1,3,6),(1,4,5),(1,4,6),(2,3,5),(2,4,5)]
Ep = [cube[k] for k in Ekeys]
webm = web(Ep)
T1, T2 = (0, 3, 5), (1, 2, 4)
h1 = intvec(rref_null([Ep[i] for i in T1], 4)[0]); h2 = intvec(rref_null([Ep[i] for i in T2], 4)[0])
mb = rref_null([h1, h2], 4); m0 = intvec(mb[0]); m1 = intvec(mb[1])
def wedge(x, y):
    return [x[i] * y[j] - x[j] * y[i] for (i, j) in combinations(range(4), 2)]
pim = wedge([Fr(c) for c in m0], [Fr(c) for c in m1])
def limdir(y, z0, w):
    # pencil R(s) = span(y, z0 + s w); R(0) = span(y, z0) = m
    yy = [P([y[i]]) for i in range(4)]; zs = [P([z0[i], w[i]]) for i in range(4)]
    p, q = yy, zs
    rows = [[bil(Qm, p, p) for Qm in webm], [bil(Qm, q, q) for Qm in webm], [bil(Qm, p, q) for Qm in webm]]
    c = pkernel(rows); B = matsum(c, webm)
    e1 = [P([x]) for x in E1]; e2 = [P([x]) for x in E2]
    T = [[bilS(B, p, e1), bilS(B, p, e2), ZERO, ZERO], [ZERO, ZERO, bilS(B, q, e1), bilS(B, q, e2)],
         [bilS(B, q, e1), bilS(B, q, e2), bilS(B, p, e1), bilS(B, p, e2)]]
    k = pkernel(T)
    w1 = [k[0] * E1[i] + k[1] * E2[i] for i in range(4)]; w2 = [k[2] * E1[i] + k[3] * E2[i] for i in range(4)]
    # Plucker derivative pi' = w1^q + p^w2, polynomial in s; divide by common s-power, evaluate at 0
    pd = [w1[i] * q[j] - w1[j] * q[i] + p[i] * w2[j] - p[j] * w2[i] for (i, j) in combinations(range(4), 2)]
    # also pi(s) itself; tangent direction = pd mod pi(s). Take limit of the 2-plane span(pi, pd) at s=0.
    pi = [p[i] * q[j] - p[j] * q[i] for (i, j) in combinations(range(4), 2)]
    # 2x2 minors of [pi; pd] (Plucker coords of the 2-plane in Lambda^2 = the tangent line of G-curve), limit at s=0
    mins = [pi[a] * pd[b] - pi[b] * pd[a] for a, b in combinations(range(6), 2)]
    g = pgcd_list(mins)
    vals = []
    for mm in mins:
        r = mm // g
        vals.append(int(r.coeffs()[0]) if r != 0 and r.degree() >= 0 else 0)
    # normalise projectively
    from math import gcd
    from functools import reduce
    gg = reduce(gcd, [abs(v) for v in vals if v], 0)
    vals = [v // gg for v in vals]
    for v in vals:
        if v != 0:
            if v < 0: vals = [-x for x in vals]
            break
    return tuple(vals)
rng = random.Random(3)
res = []
for t in range(6):
    al, be = rng.randint(1, 9), rng.randint(-9, 9); ga, de = rng.randint(-9, 9), rng.randint(1, 9)
    if al * de - be * ga == 0: continue
    y = [al * m0[i] + be * m1[i] for i in range(4)]; z0 = [ga * m0[i] + de * m1[i] for i in range(4)]
    w = [rng.randint(-20, 20) for _ in range(4)]
    d = limdir(y, z0, w)
    res.append(d)
    print('pencil', t, 'limit tangent line at m (6 Plucker minors of span(pi,pi_dot)):', d)
print('number of distinct limit directions:', len(set(res)), 'out of', len(res))
