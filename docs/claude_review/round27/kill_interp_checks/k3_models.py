"""K3: rebuild sdc7 and syl8 (L=2, T=1, n=8) from the record's definitions (round11 model8.py / bigger.py),
own code; graded Hilbert functions and (SHIFT_e) w.r.t. L_pair mod two primes."""
import sys, time
sys.path.insert(0, '.')
import numpy as np
from klib import build_box, Lpair, Lcoord
from klib2 import *
def legendre(a, q):
    a %= q
    return 0 if a == 0 else (1 if pow(a, (q - 1) // 2, q) == 1 else -1)
def canon(n, S):
    Sc = tuple(x for x in range(n) if x not in S)
    return min(tuple(S), Sc, key=lambda s: (len(s), s))
def paley_I(q):
    idx = ['inf'] + list(range(q)); M = q + 1
    H = np.zeros((M, M), dtype=int)
    for a, x in enumerate(idx):
        for b, y in enumerate(idx):
            if x == 'inf': H[a, b] = 1
            elif y == 'inf': H[a, b] = -1
            else: H[a, b] = (1 if x == y else 0) + legendre(y - x, q)
    return H
def sdc(q):
    K = paley_I(q); M = q + 1
    if not ((K + K.T) == 2 * np.eye(M, dtype=int)).all():
        K = K.T
    Dg = np.diag(K[:, 0]); K = Dg @ K @ Dg
    pats = set()
    for a in range(1, M):
        S = [x for x in range(M) if K[x, a] == -1]
        Sp = sorted(set(S) ^ {0, a})
        pats.add(canon(M, tuple(S))); pats.add(canon(M, tuple(Sp)))
    return tuple(sorted(pats))
def syl8():
    H = [[1]]
    for _ in range(3):
        H = [r + r for r in H] + [r + [-x for x in r] for r in H]
    return tuple(sorted(set(canon(8, tuple(x for x in range(8) if r[x] == 1)) for r in H[1:])))
which = sys.argv[1]; emax = int(sys.argv[2])
prof = sdc(7) if which == 'sdc7' else syl8()
A = build_box(8, [(S, 2) for S in prof])
d = 7
out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
pr(which, "profile", prof, "|A| =", len(A))
for p in (2147483629, 1000000007):
    t0 = time.time()
    Hs = []
    for e in range(0, emax + 1):
        hs, sm = shift_dims(A, Lpair(d), e, p, d)
        Hs.append(hs)
        pr("  p=%d e=%d H_A(e)=%d dim sum_{L_pair} M_e(A_l)=%d SHIFT=%s  (%.1fs)" % (p, e, hs, sm, hs == sm, time.time() - t0))
    pr("  p=%d graded H_A(0..%d) = %s" % (p, emax, Hs))
open("k3_%s_out.txt" % which, "w").write("\n".join(out) + "\n")
