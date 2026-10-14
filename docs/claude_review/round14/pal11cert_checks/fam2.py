"""rho(B_F) for 2-direction families F = (v1, v2) among the test directions with tau < 8 (exact, integer ranks).
Prop 2.6(b) of fatpoint: a product certificate has #{i : z_i' not orthogonal to span F} <= rho(B_F)."""
import numpy as np, itertools, time
from flint import fmpz_mat
from fplib import pal_patterns, box
M_, pats = pal_patterns(11); A = box(M_, pats, 2)
Aarr = np.array(A, dtype=np.int64)
dat = np.load('dirs_pal11.npz'); D = dat['D']; tau = dat['tau']
V = Aarr @ D.T   # |A| x ND
def rho_planar(B):
    B = sorted(B); n = len(B)
    rk = 0; deg = 0; rows = []
    while True:
        for a in range(deg + 1):
            rows.append([s ** (deg - a) * t ** a for (s, t) in B])
        rk_new = fmpz_mat(len(rows), n, [x for r in rows for x in r]).rank()
        if rk_new == n: return deg
        deg += 1
t0 = time.time()
out = []
cache = {}
for i in range(len(D)):
    for j in range(i + 1, len(D)):
        if tau[i] + tau[j] > 10: continue
        code = np.unique(V[:, i] * 64 + V[:, j])
        key = tuple(code.tolist())
        if key not in cache:
            B = [((c + 32) // 64, c - 64 * ((c + 32) // 64)) for c in key]
            cache[key] = rho_planar(B)
        rh = cache[key]
        if rh < 8: out.append((i, j, rh))
print('pairs (tau_i+tau_j <= 10) with rho(B_F) < 8:', len(out), 'distinct images', len(cache), '%.0fs' % (time.time() - t0))
np.save('fam2_pal11.npy', np.array(out, dtype=np.int64))
from collections import Counter
print(sorted(Counter((int(tau[i]), int(tau[j]), rh) for i, j, rh in out).items()))
