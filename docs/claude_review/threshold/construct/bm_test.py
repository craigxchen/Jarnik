"""Test 'balance monotonicity' (BM): for a finite S in Z^M, phi(r) = average over sign vectors eps with
r minus signs of width_S(eps) is nondecreasing in min(r, M-r).  Also tests the good-pair lemma (GPL)
and the covariance form.  Pure exact integer arithmetic."""
import random, itertools, sys
from fractions import Fraction as Fr
import numpy as np
random.seed(int(sys.argv[3]) if len(sys.argv) > 3 else 1)
M = int(sys.argv[1]); trials = int(sys.argv[2])
eps = np.array(list(itertools.product([1, -1], repeat=M)), dtype=np.int64)
r_of = (eps == -1).sum(axis=1)
viol_bm = 0; viol_gp = 0; worst = None
for t in range(trials):
    mode = random.choice(["pts", "few", "line"])
    if mode == "pts":
        n = random.randint(2, 12); R = random.randint(1, 6)
    elif mode == "few":
        n = random.randint(2, 4); R = random.randint(1, 20)
    else:
        n = 2; R = random.randint(1, 30)
    S = []
    for _ in range(n):
        k = [random.randint(-R, R) for _ in range(M-1)]
        S.append(k + [-sum(k)])
    A = np.array(S, dtype=np.int64)
    V = eps @ A.T                         # (2^M, n)
    W = V.max(axis=1) - V.min(axis=1)     # width along eps
    phi = {}
    for r in range(M+1):
        sel = (r_of == r)
        phi[r] = Fr(int(W[sel].sum()), int(sel.sum()))
    seq = [phi[r] for r in range(0, M//2 + 1)]
    if any(seq[i] > seq[i+1] for i in range(len(seq)-1)):
        viol_bm += 1
        if worst is None: worst = (S, seq)
    # good pair: exists i<j with E[W eps_i eps_j] <= 0
    ok = False
    for i in range(M):
        for j in range(i+1, M):
            if int((W * eps[:, i] * eps[:, j]).sum()) <= 0:
                ok = True
    if not ok: viol_gp += 1
print(f"M={M} trials={trials}: BM violations={viol_bm}, good-pair violations={viol_gp}")
if worst: print("first BM violation:", worst)
