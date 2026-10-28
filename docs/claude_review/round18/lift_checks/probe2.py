import sys, flint
from fractions import Fraction as Fr
from beta7lift import setup, emon_h, lagr4, mp
p = int(sys.argv[1])
xs = [Fr(2), Fr(3), Fr(-5), Fr(1)]
E, D = setup(xs, p)
red = lambda P: tuple(mp(c, p) for c in P)
names = {}
for lab, A, j, P in E['pts']: names[red(P)] = '.'.join(lab)
for i, P in enumerate(E['aux']): names[red(P)] = 'r%d' % (i + 1)
q = red((Fr(-1), Fr(1198, 121), Fr(-46)))
# basis of H^0(L_q): cubics vanishing at r1,r2,r3 and doubly at q -- 8-dim; compute tau_j
import beta7lift as bl
def conds(P, order):
    a, u, w = P
    v0 = [x % p for x in emon_h((a, u, w, 1), 1)]
    if order == 1: return [v0]
    T = (4*u*w % p, 4*D['kb']*a*w % p, 4*D['kc']*a*u % p)
    return [v0, bl.grad_dot(P, T, p)]
base = []
for r in D['R']: base += conds(r, 1)
base += conds(q, 1)
def kernel(rows):
    M = flint.nmod_mat(rows, p)
    K = M.nullspace()[0]
    r = M.rank(); n = 12 - r
    return [[int(K[i, c]) for i in range(12)] for c in range(n)]
H0 = kernel(base); print('h0(L_q) =', len(H0))
taus = []
for j in (1, 2, 3):
    rows = base + [c for P in D['B'][j] for c in conds(P, 2 if P == q else 1)]
    K = kernel(rows); print('j', j, 'dim', len(K)); taus.append(K[0])
# evaluate tau_j at all special points (value of cubic; zero means zero of section beyond base, except at base pts)
allpts = list(names)
for j, t in enumerate(taus):
    zs = [names[P] for P in allpts if sum(a*b for a, b in zip(t, emon_h(P + (1,), 1))) % p == 0]
    print('tau_%d zero at special pts:' % (j+1), zs)
lam = lagr4(D)
M = [[ (sum(a*b for a, b in zip(taus[j], (conds(P,2)[1] if P == q else emon_h(P + (1,), 1)))) % p) for j in range(3)] for P in D['B'][4]]
print('B4 order', [names[P] for P in D['B'][4]]); print('matrix', M)
