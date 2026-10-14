"""Exact verification of the zonotopes found by zono_dfs.py: find k0 with all 2^s subset sums in A; compute the
finite-difference measure c of w^{k0} prod (w^{v_i}-1), check it annihilates all monomials of degree < 7 on A and
that its degree-7 moment form is 7! prod <v_i,u> (checked at random integer points)."""
import itertools, random, numpy as np
from rlib import *
cases = {
 'sdc7': [(0,0,0,0,0,0,1),(0,0,0,0,0,1,-1),(0,0,1,-1,0,0,0),(0,0,1,0,0,0,0),(0,1,0,-1,0,0,0),(1,-1,0,0,0,0,0),(1,0,0,0,0,-1,0)],
 'pal7': [(0,0,0,0,0,0,1),(0,0,0,0,0,1,-1),(0,0,0,0,0,1,0),(0,0,0,0,1,0,0),(0,0,0,1,0,0,0),(0,1,0,-1,0,0,0),(0,1,0,0,-1,0,0)]}
for name, V in cases.items():
    pats = sdc(7) if name == 'sdc7' else pal(7)
    A, _ = box(8, pats, 2, 0); S = set(map(tuple, A.tolist()))
    V = [np.array(v) for v in V]
    k0s = [k for k in S if all(tuple(np.array(k) + sum(V[i] for i in I)) in S
                               for r in range(len(V)+1) for I in itertools.combinations(range(len(V)), r))]
    print(name, 'L=2 base inf: translates k0 with the whole zonotope (2^7 subset sums) in A:', len(k0s), ' e.g. k0 =', k0s[0] if k0s else None)
    k0 = np.array(k0s[0]); c = {}
    for r in range(len(V)+1):
        for I in itertools.combinations(range(len(V)), r):
            k = tuple(k0 + sum((V[i] for i in I), np.zeros(7, dtype=int)))
            c[k] = c.get(k, 0) + (-1)**(len(V) - r)
    c = {k: x for k, x in c.items() if x}
    ks = list(c); cv = np.array([c[k] for k in ks]); K = np.array(ks)
    ok = True
    for e in reduced_monos(7, 6):
        v = np.ones(len(ks), dtype=np.int64)
        for j, ej in enumerate(e):
            if ej: v *= K[:, j]**ej
        if int(np.dot(cv, v)) != 0: ok = False
    random.seed(0); okf = True
    for _ in range(10):
        u = np.array([random.randint(-9, 9) for _ in range(7)], dtype=object)
        lhs = sum(int(x) * int(np.dot(np.array(k, dtype=object), u))**7 for k, x in c.items())
        rhs = 5040
        for v in V: rhs *= int(np.dot(np.array(v, dtype=object), u))
        okf &= (lhs == rhs)
    print('   c in V(A,7) (all reduced monomials deg<=6 annihilated):', ok, ' |c|_1 =', int(np.abs(cv).sum()), ' Phi_c = 7! prod <v_i,u>:', okf)
