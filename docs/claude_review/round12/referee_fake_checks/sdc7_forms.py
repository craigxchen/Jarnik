"""sdc7, T=1, L=2, m=8: exact V(A_b,8) for every base point b (independent kernel), forms in full coordinates."""
import json, flint, numpy as np, time
from rlib import *
n = 8; q = 7; pats = sdc(q); L = 2; m = 8
out = {}
for b in range(n):
    t = time.time()
    A, others = box(n, pats, L, b)
    monos = reduced_monos(n-1, m-1)
    E = eval_matrix(A, monos)
    M = flint.fmpz_mat(E.tolist())
    X, nul = M.nullspace()
    vecs = [[int(X[i, j]) for i in range(X.nrows())] for j in range(nul)]
    assert nul == 1, nul
    c = vecs[0]
    from math import gcd
    g = 0
    for v in c: g = gcd(g, abs(v))
    c = [v//g for v in c]
    # check degree-8 moment nonzero (rho = 8): evaluate deg<=8 rank = |A|
    E8 = eval_matrix(A, reduced_monos(n-1, m))
    r8 = flint.nmod_mat((E8 % 1000000007).tolist(), 1000000007).rank()
    # full coordinates k_full (length n)
    full = np.zeros((len(A), n), dtype=np.int64); full[:, others] = A; full[:, b] = -A.sum(1)
    # even?
    idx = {tuple(k): i for i, k in enumerate(A.tolist())}
    even = all(c[idx[tuple((-np.array(k)).tolist())]] == c[i] for i, k in enumerate(A.tolist()))
    print('base', b, '|A|', len(A), 'dim V(A,8) =', nul, ' rank deg<=8 =', r8, ' |c|_inf', max(map(abs, c)), '|c|_1', sum(map(abs, c)),
          'supp', sum(1 for v in c if v), 'even', even, '%.1fs' % (time.time()-t), flush=True)
    out[b] = {'full': full.tolist(), 'c': c}
json.dump(out, open('sdc7_forms.json', 'w'))
