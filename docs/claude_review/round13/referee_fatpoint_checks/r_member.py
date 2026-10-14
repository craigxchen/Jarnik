"""Referee: membership of products of 7 heavy forms in J(Y_3(2))^perp_7, own implementation.
T_7 = { (sum_k c_k k^a)_{|a|=7} : sum_k c_k f(k) = 0 for all f of degree <= 6 }  (mod p);
P = sum_a p_a u^a is in J^perp_7 iff (p_a / multinom(7;a))_a in T_7.
Rigour: non-membership mod p => non-membership over Q PROVIDED rank_p(K_{<=7}) = rank_Q(K_{<=7}) (checked separately).
Also: exact zonotope test (k0 + sum_I lam_i g_i in Y for all I, lam_i in {2,4}) for every multiset."""
import itertools, numpy as np, pickle, sys, time
from math import factorial
from collections import Counter
from flint import nmod_mat
p = 1000000007
Y, D, T = pickle.load(open(sys.argv[1], 'rb'))
heavy, out, feas = pickle.load(open(sys.argv[2], 'rb'))
e = 7
low = [a for a in itertools.product(range(5), repeat=7) if sum(a) <= e - 1]
me = [a for a in itertools.product(range(8), repeat=7) if sum(a) == e]
Yo = Y.astype(object)
def evrows(mons):
    R = []
    for a in mons:
        col = np.ones(len(Y), dtype=object)
        for x, ax in enumerate(a):
            if ax: col = col * Yo[:, x] ** ax
        R.append([int(v) % p for v in col])
    return R
t0 = time.time()
Kl = nmod_mat(evrows(low), p)
X, nul = Kl.nullspace()
print('kernel dim (deg<=6):', nul, 'expected', len(Y) - 1597, flush=True)
V = nmod_mat(nul, len(Y), [int(X[i, k]) for k in range(nul) for i in range(len(Y))], p)
E7 = nmod_mat(evrows(me), p).transpose()
Tm = V * E7
rT = Tm.rank()
print('dim T_7 mod p =', rT, '(%.0fs)' % (time.time() - t0), flush=True)
# annihilator N of T (N t = 0): nullspace of Tm (right kernel)
Nn, nn = Tm.nullspace()
N = np.array([[int(Nn[i, k]) for i in range(Nn.nrows())] for k in range(nn)], dtype=np.int64)  # nn x |me|
print('annihilator dim', nn, flush=True)
idx = {a: i for i, a in enumerate(me)}
minv = np.array([pow(factorial(e) // int(np.prod([factorial(x) for x in a])), -1, p) for a in me], dtype=np.int64)
def prodvec(G):
    poly = {(0,) * 7: 1}
    for g in G:
        new = {}
        for ex, cf in poly.items():
            for x in range(7):
                if g[x]:
                    e2 = list(ex); e2[x] += 1; e2 = tuple(e2)
                    new[e2] = new.get(e2, 0) + cf * int(g[x])
        poly = {k: v for k, v in new.items() if v}
    vec = np.zeros(len(me), dtype=np.int64)
    for ex, cf in poly.items(): vec[idx[ex]] = (cf % p) * int(minv[idx[ex]]) % p
    return vec
Yset = set(map(tuple, Y.tolist()))
subs = np.array(list(itertools.product((0, 1), repeat=7)))
def zonofits(G):
    for lams in itertools.product((2, 4), repeat=7):
        Z = np.array([l * g for l, g in zip(lams, G)])
        S = subs @ Z
        lo = -2 - S.min(0); hi = 2 - S.max(0)
        if (lo > hi).any(): continue
        for k0 in itertools.product(*[range(l, h + 1) for l, h in zip(lo, hi)]):
            if all(tuple(r) in Yset for r in (S + np.array(k0)).tolist()): return True
    return False
stats = Counter(); t0 = time.time(); bad = []
for n_, ms in enumerate(out):
    G = [heavy[i] for i in ms]
    v = prodvec(G)
    mem = not ((N @ v) % p).any() if False else not (np.array([(int(r @ v) % p) for r in []]).any())
    # exact-safe modular dot products: split to avoid int64 overflow
    acc = np.zeros(nn, dtype=np.int64)
    for lo_ in range(0, len(me), 64):
        acc = (acc + (N[:, lo_:lo_ + 64] * v[lo_:lo_ + 64]).sum(1) % p) % p
    mem = not acc.any()
    ax = int(np.max(np.sum([g != 0 for g in G], 0)))
    zf = zonofits(G) if ax <= 2 or mem else None
    stats[(ax, mem, zf)] += 1
    if mem and not zf: bad.append(ms)
print('(maxaxis, member mod p, zonotope fits) -> count   (%.0fs)' % (time.time() - t0))
for k, v in sorted(stats.items(), key=str): print('  ', k, v)
print('members that are not zonotopes:', len(bad), bad[:5])
