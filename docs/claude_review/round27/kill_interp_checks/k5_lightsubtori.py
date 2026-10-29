"""K5: in M8 (= syl8, L = 2), test whether the codim-1 subtorus K_v = {Y^v = 1} lies in Bs(V(A,m)) for the
non-pair binomial characters v found by round 11 (had8_T1_L2_m8_cone_out.txt).
K_v subset Bs(V(A,m))  <=>  every function on A that is constant on v-fibres (k ~ k + v) is in Poly_{<m}|_A.
Test: rank[Poly_{<m}|_A ; fibre indicators] == rank Poly_{<m}|_A   (mod p, two primes; and exact for the decisive m).
A generic point of K_v has Y^l != 1 for every l in L_pair (checked: v not a multiple of any l in L_pair and the
intersection K_v cap {Y^l = 1} is a proper subtorus of K_v), so K_v subset Bs(V(A,m)) forces (SHIFT_e) w.r.t.
L_pair to fail for some e < m (interp Thm 3.1(1), contrapositive)."""
import sys, itertools
sys.path.insert(0, '.')
import numpy as np, flint
from klib import build_box, monos_lt, Lpair
from klib2 import evalmat, nmod
def canon(n, S):
    Sc = tuple(x for x in range(n) if x not in S)
    return min(tuple(S), Sc, key=lambda s: (len(s), s))
def syl8():
    H = [[1]]
    for _ in range(3):
        H = [r + r for r in H] + [r + [-x for x in r] for r in H]
    return tuple(sorted(set(canon(8, tuple(x for x in range(8) if r[x] == 1)) for r in H[1:])))
A = build_box(8, [(S, 2) for S in syl8()])
d = 7
chars = {'w0 - w1*w6': {0: 1, 1: -1, 6: -1}, 'w0 - w2*w5': {0: 1, 2: -1, 5: -1}, 'w0*w3 - w4': {0: 1, 3: 1, 4: -1},
         'w0*w5 - w2': {0: 1, 5: 1, 2: -1}, 'w0*w6 - w1': {0: 1, 6: 1, 1: -1}, 'w1*w3 - w5': {1: 1, 3: 1, 5: -1},
         'w2*w3 - w6': {2: 1, 3: 1, 6: -1}}
out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
Aset = set(A); idx = {k: i for i, k in enumerate(A)}
pr("|A| =", len(A))
for name, vd in chars.items():
    v = tuple(vd.get(x, 0) for x in range(d))
    # fibres of k ~ k + v
    seen = set(); fibres = []
    for k in A:
        if k in seen: continue
        # walk to the start of the fibre
        s = k
        while tuple(a - b for a, b in zip(s, v)) in Aset: s = tuple(a - b for a, b in zip(s, v))
        fib = []
        t = s
        while t in Aset:
            fib.append(t); seen.add(t); t = tuple(a + b for a, b in zip(t, v))
        fibres.append(fib)
    nontriv = sum(1 for f in fibres if len(f) > 1)
    res = []
    for p in (2147483629, 1000000007):
        for m in (5, 6, 7, 8):
            P = evalmat(A, monos_lt(d, m), p)
            r0 = nmod(P, p).rank()
            Ind = np.zeros((len(fibres), len(A)), dtype=np.int64)
            for i, f in enumerate(fibres):
                for k in f: Ind[i, idx[k]] = 1
            r1 = nmod(np.vstack([P, Ind]), p).rank()
            res.append((p % 1000, m, r0, r1, r0 == r1))
    pr("%-12s v=%s #fibres=%d (nontrivial %d): (p mod 1000, m, rank Poly_<m, rank with fibre fns, K_v in Bs) %s" % (
        name, v, len(fibres), nontriv, res))
open("k5_lightsubtori_out.txt", "w").write("\n".join(out) + "\n")
