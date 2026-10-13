# Independent recomputation of Proposition 5.4: sum over the S_7-orbit of classes of comp(sigma f_x).
import sys, time; sys.path.insert(0, '.')
from reflib import *
from fractions import Fraction as Fr
x = tuple(int(a) for a in sys.argv[1].split(',')) if len(sys.argv) > 1 else (6, -6, -3, -5, 2, 0, 1)
cube = {(i,j,k): meet3(i,j,k) for i in (1,2) for j in (3,4) for k in (5,6)}
Ekeys = [(1,3,5),(1,3,6),(1,4,5),(1,4,6),(2,3,5),(2,4,5)]
webm = web([cube[k] for k in Ekeys])
V, data = member_through(x, webm)
cls0, good0, _ = class_and_goodness(V)
assert good0
def psi_dot(i, cls):
    j, k = [a for a in range(7) if a != i][:2]
    tot = 0
    for S, d in cls.items():
        # D_S = D_{S^c}; T ranges over subsets with i in T, j,k not in T
        Sc = frozenset(range(7)) - S
        for T in (S, Sc):
            if i in T and j not in T and k not in T:
                tot += d
    return tot
def lam_forget(i, cls):
    return -psi_dot(i, cls) + sum(cls.get(frozenset((i, j)), 0) for j in range(7) if j != i)
seen = {}
t0 = time.time()
nbad = 0
for sg in permutations(range(7)):
    W = [None] * 7
    for a in range(7):
        W[sg[a]] = V[a]
    # class of relabelled curve: D_S.(sigma f) = D_{sigma^{-1} S}.f
    cls = {frozenset(sg[a] for a in S): d for S, d in cls0.items()}
    key = class_vector(cls)
    if key in seen:
        continue
    Br, dBr, _ = FA_brackets(W, webm)
    c, nk = comp_count(Br, dBr)
    # forgetful validation on the same curve (Lemma 4.1 exact vs Lemma 5.3 count)
    fv = []
    for i in range(7):
        Bf, dBf = forget_brackets(W, i, (1, 0))
        cf, _ = comp_count(Bf, dBf)
        fv.append((cf, lam_forget(i, cls)))
        if cf != lam_forget(i, cls):
            nbad += 1
    seen[key] = (c, nk, sg, fv)
print('x =', x, ' distinct classes:', len(seen), ' time %.0fs' % (time.time() - t0))
vals = [v[0] for v in seen.values()]
S = sum(vals)
print('sum comp =', S, ' avg =', Fr(S, len(vals)), ' min', min(vals), ' max', max(vals))
from collections import Counter
print('value distribution:', sorted(Counter(vals).items()))
print('forgetful validation mismatches (Lemma 4.1 vs count), over 7*#classes:', nbad)
print('min #nonzero W_k:', min(v[1] for v in seen.values()))
import pickle
pickle.dump({k: (v[0], v[2]) for k, v in seen.items()}, open('t3_vals_%s.pkl' % ('_'.join(map(str, x))), 'wb'))
