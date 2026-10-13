# (1) parameter space of V_A: double cover of P(Lambda) branched along det=0; irreducible iff det is not a square (times const).
# (2) lambda_{f_i}(tau_* gamma_A) value set via Lemma 4.1.
import sys; sys.path.insert(0, '.')
from reflib import *
import sympy as sp
cube = {(i,j,k): meet3(i,j,k) for i in (1,2) for j in (3,4) for k in (5,6)}
Ekeys = [(1,3,5),(1,3,6),(1,4,5),(1,4,6),(2,3,5),(2,4,5)]
webm = web([cube[k] for k in Ekeys])
c = sp.symbols('c0:4')
M = sp.zeros(4, 4)
for k in range(4):
    M += c[k] * sp.Matrix(webm[k])
d = sp.factor(M.det())
print('det of web quadric factors as:', sp.factor_list(d))
# (2)
GA = {('O','a1','a2'): 2, ('O','b1','b2'): 2, ('O','c1','c2'): 2, ('a1','b1','c1'): 1, ('a1','b1','c2'): 1,
      ('a1','b2','c1'): 1, ('a1','b2','c2'): 1, ('a2','b1','c1'): 1, ('a2','b2','c1'): 1, ('a2','b1'): 1,
      ('a2','b2'): 1, ('a2','c2'): 2, ('b1','c2'): 1, ('b2','c2'): 1}
cls0 = {frozenset(IDX[l] for l in k): v for k, v in GA.items()}
def lamf(i, cls):
    j, k = [a for a in range(7) if a != i][:2]
    psi = 0
    for S, d_ in cls.items():
        for T in (S, frozenset(range(7)) - S):
            if i in T and j not in T and k not in T: psi += d_
    return -psi + sum(cls.get(frozenset((i, b)), 0) for b in range(7) if b != i)
vals = set()
for sg in permutations(range(7)):
    cls = {frozenset(sg[a] for a in S): d_ for S, d_ in cls0.items()}
    vals.add(lamf(0, cls))
print('lambda_{f_i}(tau gamma_A) values:', sorted(vals))
# beta check of Lemma 4.1
bcls = {frozenset(S): 1 for r in (2, 3) for S in combinations(range(7), r)}
print('lambda_{f_i}(beta):', [lamf(i, bcls) for i in range(7)])
