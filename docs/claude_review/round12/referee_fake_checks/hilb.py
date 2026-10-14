import sys, time, flint, numpy as np
from rlib import *
P1, P2 = 1000000007, 998244353
def hilbert(A, d, maxdeg, p):
    monos = reduced_monos(d, maxdeg)
    E = eval_matrix(A, monos) % p
    ranks = []
    for j in range(maxdeg+1):
        rows = [i for i, e in enumerate(monos) if sum(e) <= j]
        M = flint.nmod_mat(E[rows].tolist(), p)
        ranks.append(M.rank())
    return [ranks[0]] + [ranks[j]-ranks[j-1] for j in range(1, len(ranks))]
models = {'sdc7': (8, sdc(7)), 'pal7': (8, pal(7)), 'pal7core': (8, pal(7, core_only=True)), 'syl8': (8, sylvester(3))}
for name, (n, pats) in models.items():
    for L in (2, 3):
        for base in (0, n-1):
            t = time.time(); A, _ = box(n, pats, L, base)
            H1 = hilbert(A, n-1, 12, P1); H2 = hilbert(A, n-1, 12, P2)
            print(name, '#pats', len(pats), 'L', L, 'base', base, '|A|', len(A), 'H', H1, 'same2p', H1 == H2, 'sum', sum(H1), '%.1fs' % (time.time()-t), flush=True)
