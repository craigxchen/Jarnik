"""Generator pools for csearch3 (pal11, L=2, m=8).  Same file format as mkpool.py.
usage: python3 mkpool2.py MODE out.txt
MODE  'bal:S:E'  balanced primitive z in Z^12 (row 0 = base), support 2..S, |z_x| <= E for all 12 rows;
      'tern'     all primitive z' in {-1,0,1}^11 (z_0 = -sum z' unrestricted);
      'tern2'    all z' in {-2..2}^11 with support <= 6 (primitive);
Every generator is sign-normalised (first nonzero entry of z positive); duplicates removed."""
import sys, itertools
from math import gcd
import numpy as np
mode, out = sys.argv[1], sys.argv[2]
dat = np.load('dirs_pal11.npz'); D = dat['D']; tau = dat['tau']
n = 12
def norm(z):
    f = next(x for x in z if x)
    return tuple(z) if f > 0 else tuple(-x for x in z)
def prim(z):
    g = 0
    for v in z: g = gcd(g, abs(v))
    return g == 1
S = set()
if mode.startswith('bal'):
    _, ms, me = mode.split(':'); ms, me = int(ms), int(me)
    for s in range(2, ms + 1):
        for supp in itertools.combinations(range(n), s):
            for vals in itertools.product([e for e in range(-me, me + 1) if e], repeat=s):
                if sum(vals) != 0 or vals[0] < 0 or not prim(vals): continue
                z = [0]*n
                for x, v in zip(supp, vals): z[x] = v
                S.add(tuple(z))
elif mode == 'tern':
    for zp in itertools.product((-1, 0, 1), repeat=11):
        if any(zp):
            z = (-sum(zp),) + zp
            S.add(norm(z))
elif mode == 'tern2':
    for s in range(1, 7):
        for supp in itertools.combinations(range(11), s):
            for vals in itertools.product((-2, -1, 1, 2), repeat=s):
                if not prim(vals): continue
                zp = [0]*11
                for x, v in zip(supp, vals): zp[x] = v
                S.add(norm((-sum(zp),) + tuple(zp)))
from fplib import pal_patterns
_M, _pats = pal_patterns(11)
_core = []
for _P in _pats:
    if len(_P) == 6:
        _v = [0]*11
        for _x in _P: _v[_x-1] = 1
        _core.append(_v)
_C = np.array(_core, dtype=np.int64).T
def tcore(z): return int(((np.array(z[1:]) @ _C) != 0).sum())
# order: number of touched core columns first (exceptional generators, t_core <= 5, come first), then support
pool = sorted(S, key=lambda z: (tcore(z), sum(1 for x in z if x), sum(abs(x) for x in z), z))
index = {z: i for i, z in enumerate(pool)}
def act(a, z):
    w = [0]*n; w[0] = z[0]
    for i in range(11): w[1 + (a*i) % 11] = z[1+i]
    return norm(w)
rep = [1 if i == min(index[act(a, z)] for a in (1, 3, 4, 5, 9)) else 0 for i, z in enumerate(pool)]
P = np.array(pool, dtype=np.int64)[:, 1:]
with open(out, 'w') as f:
    f.write(f'{len(pool)} {len(D)}\n')
    f.write(' '.join(str(int(t)) for t in tau) + '\n')
    for s0 in range(0, len(pool), 5000):
        NZ = (P[s0:s0+5000] @ D.T) != 0
        for ii in range(NZ.shape[0]):
            i = s0 + ii
            idx = np.nonzero(NZ[ii])[0]
            f.write(f'{rep[i]} ' + ' '.join(map(str, pool[i])) + f' {len(idx)} ' + ' '.join(map(str, idx)) + '\n')
print(mode, 'pool', len(pool), 'orbit reps', sum(rep))
