"""Own construction of the ordinary Paley double-core rows; projected affine ranks for
anchor triples (III.4 COMPUTED claim) and a numerical check of the edge rank inequality
s*D/2 - 1 <= 2*rank (exact III.4 statement) plus the operator norm bound ||C||^2 <= 2."""
from itertools import combinations
import numpy as np
import sympy as sp

def rows_paley(q):
    sq = {x*x % q for x in range(1,q)}
    chi = lambda x: 0 if x % q == 0 else (1 if x % q in sq else -1)
    cols = []
    for a in range(q):
        for eps in (1,0):
            cols.append([eps if x == a else (1 if chi(a-x) == 1 else 0) for x in range(q)])
    return [[c[x] for c in cols] for x in range(q)]  # q rows, 2q binary columns

def affrank(vecs):
    base = vecs[0]
    M = sp.Matrix([[a-b for a,b in zip(v,base)] for v in vecs[1:]])
    return M.rank()

for q in (7,11,19):
    R = rows_paley(q); r = 2*q
    # pair separation in equal-weight model: distance > r/2 ?
    sep = min(sum(abs(a-b) for a,b in zip(R[x],R[y])) for x,y in combinations(range(q),2))
    ranks = {}
    for i,j,k in combinations(range(q),3):
        cols = [c for c in range(r) if len({R[i][c],R[j][c],R[k][c]}) == 2]
        others = [x for x in range(q) if x not in (i,j,k)]
        rk = affrank([[R[x][c] for c in cols] for x in others])
        ranks[rk] = ranks.get(rk,0)+1
    print('q',q,'min pair distance',sep,'of',r,'columns; triple projected affine ranks',ranks)
    # edge version
    worst = None
    for i,j in combinations(range(q),2):
        cols = [c for c in range(r) if R[i][c] != R[j][c]]
        D = len(cols)/r
        others = [x for x in range(q) if x not in (i,j)]
        s = len(others)
        # c_x(p) = (2 z - L) sqrt(log p/(W L)); L=1, equal weights -> (2z-1)/sqrt(r)
        C = np.array([[(2*R[x][c]-1)/np.sqrt(r) for c in cols] for x in others])
        Cc = C - C.mean(axis=0)
        tot = (Cc**2).sum()
        rk = affrank([[R[x][c] for c in cols] for x in others])
        opn = np.linalg.norm(C, 2)**2
        assert tot >= s*D/2 - 1 - 1e-9, (tot, s*D/2-1)
        assert opn <= 2 + 1e-9
        assert tot <= 2*rk + 1e-9
        slack = 2*rk - (s*D/2-1)
        worst = slack if worst is None else min(worst, slack)
    print('   edge version: all pairs satisfy sD/2-1 <= total centered norm <= 2*rank, ||C||^2<=2; min slack 2rank-(sD/2-1) =', round(worst,3))
