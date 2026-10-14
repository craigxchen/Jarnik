"""sdc11, T=1, L=2: dim V(A_b, 8) mod p for base b = 0 (inf) and b = 1, and for the common core C (all |k_x|<=1)."""
import sys, time, flint, numpy as np
from rlib import *
n = 12; pats = sdc(11); L = 2; p = 1000000007
which = sys.argv[1]
A0, o0 = box(n, pats, L, 0)
full0 = np.zeros((len(A0), n), dtype=np.int64); full0[:, o0] = A0; full0[:, 0] = -A0.sum(1)
core = full0[np.abs(full0).max(1) <= 1]
print('|A_0| =', len(A0), ' |core C| =', len(core), flush=True)
if which == 'core':
    S = core[:, 1:]          # coordinates relative to base 0 (k_0 determined)
elif which == 'b0':
    S = A0
else:
    b = int(which); Ab, ob = box(n, pats, L, b); S = Ab
    print('|A_%d| =' % b, len(S), flush=True)
t = time.time()
monos = reduced_monos(n-1, 7)
E = eval_matrix(S, monos) % p
M = flint.nmod_mat(E.tolist(), p); del E
r = M.rank()
print(which, ': |S| =', len(S), ' rank E_{<=7} mod p =', r, ' => dim V(S,8) mod p =', len(S) - r, ' (%.0fs)' % (time.time()-t), flush=True)
