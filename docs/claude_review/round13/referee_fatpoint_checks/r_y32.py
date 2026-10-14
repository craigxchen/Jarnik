"""Referee: independent construction of Y_3(2), widths tau_Y(v), and enumeration of heavy multisets.
Y_3(2) := Lambda_U ∩ [-2,2]^7, Lambda_U = { y : c_zeta = -(1/4) sum_{a!=0} y_a (-1)^{a.zeta} integral for all zeta }.
(derived from y = M c, c in Z^8_0, M_{a,zeta} = [a.zeta=1]).
"""
import itertools, numpy as np, pickle, sys
axes = list(range(1, 8))
def chi(a, z): return -1 if bin(a & z).count('1') % 2 else 1
CH = np.array([[chi(a, z) for a in axes] for z in range(8)])  # 8 x 7
pts = np.array(list(itertools.product(range(-2, 3), repeat=7)))
S = pts @ CH.T  # sum_a y_a chi_z(a)
ok = (S % 4 == 0).all(axis=1)
Y = pts[ok]
print('|Y| =', len(Y))
# cross-check: y = M c with c integral and sum c = 0
c = -(S[ok]) // 4
M = np.array([[1 if bin(a & z).count('1') % 2 else 0 for z in range(8)] for a in axes])
assert (c @ M.T == Y).all() and (c.sum(axis=1) == 0).all()
# widths over directions v in {-2..2}^7 (primitive, up to sign)
dirs = []
for v in itertools.product(range(-2, 3), repeat=7):
    v = np.array(v)
    if not v.any(): continue
    if np.gcd.reduce(np.abs(v[v != 0])) != 1: continue
    fnz = v[np.nonzero(v)[0][0]]
    if fnz < 0: continue
    dirs.append(v)
dirs = np.array(dirs)
vals = Y @ dirs.T
tau = np.array([len(np.unique(vals[:, i])) - 1 for i in range(len(dirs))])
cheap = tau < 8
print('directions tested', len(dirs), ' with tau<=6:', (tau <= 6).sum(), ' tau==7:', (tau == 7).sum())
from collections import Counter
print('tau histogram (<=7):', sorted(Counter(tau[tau <= 7]).items()))
for i in np.nonzero(tau <= 6)[0][:60]:
    print('  v=', dirs[i].tolist(), 'tau=', tau[i], 'values', np.unique(vals[:, i]).tolist())
pickle.dump((Y, dirs[tau <= 7], tau[tau <= 7]), open(sys.argv[1], 'wb'))
