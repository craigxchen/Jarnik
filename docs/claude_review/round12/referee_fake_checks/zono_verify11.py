import itertools, random, numpy as np
from rlib import *
V = [(1,0,0,0,0,0,0,0,0,0,0),(0,1,0,0,0,0,0,0,0,0,0),(1,-1,0,0,0,0,0,0,0,0,0),(0,0,1,0,0,-1,0,0,0,0,0),
     (0,0,1,0,0,0,0,0,0,0,-1),(0,0,0,0,0,1,0,-1,0,0,0),(0,0,0,0,0,0,0,1,0,0,-1)]
A, _ = box(12, sdc(11), 2, 0); S = set(map(tuple, A.tolist())); V = [np.array(v) for v in V]
subs = [sum((V[i] for i in I), np.zeros(11, dtype=int)) for r in range(8) for I in itertools.combinations(range(7), r)]
k0s = [k for k in S if all(tuple(np.array(k) + s) in S for s in subs)]
print('sdc11 L=2 base inf: translates with all 128 subset sums in A:', len(k0s), 'e.g.', k0s[0])
k0 = np.array(k0s[0]); c = {}
for r in range(8):
    for I in itertools.combinations(range(7), r):
        k = tuple(k0 + sum((V[i] for i in I), np.zeros(11, dtype=int))); c[k] = c.get(k, 0) + (-1)**(7 - r)
c = {k: x for k, x in c.items() if x}; ks = list(c); cv = np.array([c[k] for k in ks]); K = np.array(ks)
ok = all(int(np.dot(cv, np.prod([K[:, j]**e[j] for j in range(11)], axis=0))) == 0 for e in reduced_monos(11, 6))
random.seed(1); okf = True
for _ in range(10):
    u = np.array([random.randint(-9, 9) for _ in range(11)], dtype=object)
    lhs = sum(int(x) * int(np.dot(np.array(k, dtype=object), u))**7 for k, x in c.items())
    rhs = 5040
    for v in V: rhs *= int(np.dot(np.array(v, dtype=object), u))
    okf &= lhs == rhs
print('c in V(A,7):', ok, ' |c|_1 =', int(np.abs(cv).sum()), ' Phi_c = 7! prod <v_i,u>:', okf)
