# Locus of lines R with rank(Lambda|_R) <= 2 (web quadric through R not unique): sample a 2-plane section.
import sys, random; sys.path.insert(0, '.')
from reflib import *
import sympy as sp
cube = {(i,j,k): meet3(i,j,k) for i in (1,2) for j in (3,4) for k in (5,6)}
Ekeys = [(1,3,5),(1,3,6),(1,4,5),(1,4,6),(2,3,5),(2,4,5)]
webm = web([cube[k] for k in Ekeys])
u, v = sp.symbols('u v')
rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 3)
co = [[rng.randint(-5, 5) for _ in range(3)] for _ in range(4)]
a, b, c, d = [co[i][0] * u + co[i][1] * v + co[i][2] for i in range(4)]
p = [1, 0, a, b]; q = [0, 1, c, d]
def B(S, x, y):
    return sp.expand(sum(S[i][j] * x[i] * y[j] for i in range(4) for j in range(4)))
M = sp.Matrix([[B(S, p, p), B(S, p, q), B(S, q, q)] for S in webm])
minors = [sp.expand(M.extract([r for r in range(4) if r != k], [0, 1, 2]).det()) for k in range(4)]
G = sp.groebner(minors, u, v, order='lex')
print('groebner basis lengths', len(G.exprs), [sp.Poly(g, u, v).total_degree() for g in G.exprs])
last = sp.Poly(G.exprs[-1], v)
print('univariate in v: degree', last.degree(), ' factors:', [(sp.Poly(f, v).degree(), e) for f, e in sp.factor_list(G.exprs[-1])[1]])
sols = sp.solve(G.exprs, [u, v], dict=True)
print('#solutions', len(sols))
for s_ in sols:
    uu, vv = s_[u], s_[v]
    P_ = [sp.nsimplify(sp.sympify(x).subs({u: uu, v: vv})) for x in p]; Q_ = [sp.nsimplify(sp.sympify(x).subs({u: uu, v: vv})) for x in q]
    # classify: through which cube points; contained in which planes; which planes meet R at same point
    thr = [k for k in cube if sp.Matrix([P_, Q_, cube[k]]).rank() == 2]
    inpl = [LAB[i] for i in range(7) if sp.simplify(sum(PL[i][j] * P_[j] for j in range(4))) == 0 and sp.simplify(sum(PL[i][j] * Q_[j] for j in range(4))) == 0]
    # positions x_a = (l_a(q), -l_a(p))
    pos = []
    for i in range(7):
        lp = sp.simplify(sum(PL[i][j] * P_[j] for j in range(4))); lq = sp.simplify(sum(PL[i][j] * Q_[j] for j in range(4)))
        pos.append((lq, -lp))
    coll = [(LAB[i], LAB[j]) for i, j in combinations(range(7), 2) if sp.simplify(pos[i][0] * pos[j][1] - pos[i][1] * pos[j][0]) == 0]
    print('sol: thr cube pts', [''.join(LAB[i] for i in k) for k in thr], ' in planes', inpl, ' collisions', coll,
          ' rational:', all(x.is_rational for x in P_ + Q_))
