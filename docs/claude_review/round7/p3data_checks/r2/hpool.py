"""Semi-continuous MILP (qsolve.sc_milp) over a structured pool: columns of K random row-permutations
of the Paley Hadamard design (M-1 balanced columns each).  Upper bounds for W_sc only.
Usage: python3 hpool.py M t Xf K time_limit [extra_jsonl]"""
import sys, math, time, json
import numpy as np
sys.path.insert(0, '..')
from gauss7 import paley_normalised
from qcore import primes_1mod4_above, s_to_mask
from qsolve import sc_milp

def pool_masks(M, K, rng, extra=()):
    H = paley_normalised(M - 1)[:, 1:]
    out = set(extra)
    for k in range(K):
        perm = np.arange(M) if k == 0 else rng.permutation(M)
        Hp = H[perm]
        for j in range(M - 1):
            out.add(s_to_mask(Hp[:, j]))
    return sorted(out)

if __name__ == '__main__':
    M, t, Xf, K, tl = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), float(sys.argv[5])
    D = np.load('../data/D_construct_M%d_t%d_u1_X%d.npy' % (M, t, Xf))
    X = M * Xf; w0 = math.log(primes_1mod4_above(X, 1)[0])
    rng = np.random.default_rng(t)
    masks = pool_masks(M, K, rng)
    t0 = time.time()
    out = sc_milp(D, masks, w0, time_limit=tl, gap=1e-4)
    lam = out['lam']
    print('M=%d t=%d X=%d pool %d: W %s (cols %s) mip lb %s  %.0fs' % (M, t, X, len(masks), out['W'],
          None if lam is None else int((lam > 0).sum()), out['lb'], time.time() - t0), flush=True)
