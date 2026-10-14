import time, flint, numpy as np
from rlib import *
p = 1000000007
A, _ = box(12, pal(11), 2, 0)
t = time.time()
for deg in (7, 8):
    E = eval_matrix(A, reduced_monos(11, deg)) % p
    r = flint.nmod_mat(E.tolist(), p).rank(); del E
    print('pal11 L=2 base inf |A| =', len(A), ' rank Poly_{<=%d}|A mod p =' % deg, r, ' => dim V(A,%d) mod p = %d' % (deg+1, len(A)-r), ' (%.0fs)' % (time.time()-t), flush=True)
