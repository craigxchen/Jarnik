"""Exact checks of the bookkeeping in VI.5 on random affine pencils G_e = a_e B + c_e:
 (1) t_ij = [P_i,P_j]/|G_ij|^2 is an integer and equals Im(kappa * conj(A(B)) * D(B)),
     kappa = conj(prod a over i-star minus ij) * prod a over j-star minus ij, A,D the monic root cubics;
 (2) Im(D(B)/(A-D)(B)) = t/(kappa |A-D|^2) when kappa is real (here forced by a_e real);
 (3) distribution of m = deg(A-D) over the ten pairs, and the size of Im Q(B) for generic B
     (it is NOT zero generically: the descent only applies when |t| <= P^(nu)).
"""
import random, sys
import sympy as sp
from itertools import combinations
sys.path.insert(0, '.')
from gauss import EDGES, cofactors
I = sp.I
u = sp.Symbol('u')
random.seed(3)
for trial in range(4):
    a = [random.choice([1, 2, 3, -1, -2]) for _ in range(10)]          # real a_e (aligned phases)
    c = [random.randint(-9, 9) + I*random.choice([-5,-4,-3,-2,-1,1,2,3,4,5]) for _ in range(10)]
    B = random.randint(10**6, 2*10**6) + I*random.randint(-10**6, 10**6)
    G = [sp.expand(ae*B + ce) for ae, ce in zip(a, c)]
    Gi = [(int(sp.re(g)), int(sp.im(g))) for g in G]
    t = cofactors(Gi)
    roots = [-ce/ae for ae, ce in zip(a, c)]
    ms = []
    for (i, j) in combinations(range(5), 2):
        Ai = [k for k, e in enumerate(EDGES) if i in e and e != (i, j)]
        Dj = [k for k, e in enumerate(EDGES) if j in e and e != (i, j)]
        A = sp.expand(sp.prod([u - roots[k] for k in Ai]))
        D = sp.expand(sp.prod([u - roots[k] for k in Dj]))
        kappa = sp.prod([a[k] for k in Ai]) * sp.prod([a[k] for k in Dj])   # real here
        AB = sp.expand(A.subs(u, B)); DB = sp.expand(D.subs(u, B))
        val = sp.im(sp.expand(kappa*sp.conjugate(AB)*DB))
        assert sp.simplify(val - t[(i, j)]) == 0, (val, t[(i, j)])
        diff = sp.Poly(sp.expand(A - D), u)
        m = diff.degree(); ms.append(m)
        if m >= 1:
            Q, R = sp.div(sp.Poly(D, u), diff)
            lhs = sp.im(sp.nsimplify(sp.expand(DB / diff.as_expr().subs(u, B))))
            rhs = sp.Rational(t[(i, j)]) / (kappa * sp.Abs(sp.expand(diff.as_expr().subs(u, B)))**2)
            assert sp.simplify(lhs - rhs) == 0
            imQ = sp.im(sp.expand(Q.as_expr().subs(u, B)))
    print(f'trial {trial}: identities (1),(2) exact for all ten pairs; m-distribution {sorted(ms)}; '
          f'|t| range [{min(abs(v) for v in t.values()):.2e}, {max(abs(v) for v in t.values()):.2e}], last |Im Q(B)| = {float(abs(imQ)):.3e}')
