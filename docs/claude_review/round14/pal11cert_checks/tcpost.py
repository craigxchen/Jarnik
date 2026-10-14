"""(1) EXACT: coordinate subspaces span(e_S) contained in TC(A,8).  For S a set of coordinates of k', B_S = pi_S(A).
    If rank_p K_{<=7}(B_S) = |B_S| then rank_Q = |B_S|, i.e. rho(B_S) <= 7, J(B_S)_8 = all forms, and by fatpoint
    Prop 2.2(b) + Lemma 2.1, span(e_S) = Bs(V_{F_S}(8)) is contained in TC(A,8).  All S with |S| <= 8 are tested.
(2) EVIDENCE: the numerical zeros of tcnum.py (pickle): support, Jacobian rank, containing coordinate subspace."""
import pickle, itertools, sys
import numpy as np
from flint import nmod_mat
from fplib import pal_patterns, box, reduced_monomials
p = 2147483629
M_, pats = pal_patterns(11); A = box(M_, pats, 2); d = 11
Aarr = np.array(A, dtype=np.int64)
monos = {s: [a for e in range(8) for a in reduced_monomials(s, e)] for s in range(1, 9)}
def rho_lt8(S):
    B = np.unique(Aarr[:, list(S)], axis=0)
    mons = monos[len(S)]
    if len(mons) < len(B): return False, len(B)
    E = np.ones((len(mons), len(B)), dtype=np.int64)
    for r, a in enumerate(mons):
        for x, ax in enumerate(a):
            if ax == 1: E[r] *= B[:, x]
            elif ax == 2: E[r] *= B[:, x] ** 2
    rk = nmod_mat(E.shape[0], E.shape[1], [int(v) % p for v in E.flatten()], p).rank()
    return rk == len(B), len(B)
inTC = {}
maxsets = []
for s in range(1, 9):
    cnt = 0
    for S in itertools.combinations(range(d), s):
        ok, nb = rho_lt8(S)
        if ok:
            inTC[S] = True; cnt += 1
    print('|S| =', s, ': coordinate subspaces span(e_S) certified inside TC(A,8):', cnt, flush=True)
# maximal ones
allS = sorted(inTC, key=len, reverse=True)
maximal = [S for S in allS if not any(set(S) < set(T) for T in allS if len(T) > len(S))]
from collections import Counter
print('maximal certified coordinate subspaces by dimension:', sorted(Counter(len(S) for S in maximal).items()))
rel = {tuple(x - 1 for x in P): P for P in pats}
for S in maximal:
    if len(S) >= 5:
        tag = 'pattern' if S in rel else ''
        print('   maximal', S, tag)
pickle.dump(maximal, open('tc_coord_maximal.pkl', 'wb'))
# (2) numerical zeros
zeros = pickle.load(open(sys.argv[1], 'rb'))
stats = Counter()
for F, u in zeros:
    u = u / np.abs(u).max()
    S = tuple(np.nonzero(np.abs(u) > 1e-7)[0].tolist())
    inside = any(set(S) <= set(T) for T in maximal)
    stats[(len(S), inside, F < 1e-16)] += 1
print('numerical zeros: (support size, support inside a certified coordinate subspace of TC, F<1e-16) -> count')
for k, v in sorted(stats.items()): print('   ', k, v)
