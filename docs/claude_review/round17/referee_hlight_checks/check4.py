import numpy as np, itertools
from rlib import *
out = []
def log(s): print(s, flush=True); out.append(s)
# brute-force 4-term F=0 minimum WITHOUT symmetry reduction, from the pattern definition
for q in [19, 23, 31]:
    P = patterns(q); n = q+1; best = 10**9
    for quad in itertools.combinations(range(1, n), 4):
        for plus in itertools.combinations(quad, 2):
            lam = np.zeros(n, dtype=np.int64); lam[list(quad)] = -1; lam[list(plus)] = 1
            best = min(best, sigma_def(P, lam)//2)
    log(f"q={q}: brute-force (no symmetry) 4-term F=0 min sigma = {best}")
# dual pairs from the definition (small q) and via exact convolution (large q)
def dual(q):
    ch = chi_table(q); f = np.zeros(q, dtype=np.int64)
    for x in range(2, q): f[x] = (ch[x] - ch[(x-1) % q]) // 2
    f[1] = 1; return f
for q in [7, 11, 19, 23, 31, 43, 83]:
    P = patterns(q); f = dual(q); lam = np.concatenate([[0], f])
    log(f"q={q}: dual pair sum={lam.sum()} l1={np.abs(lam).sum()} ((q+1)/2={(q+1)//2}) sigma(def)={sigma_def(P, lam)/2} ((q-1)/2={(q-1)/2})")
for q in [3251, 3259, 5479]:
    ch = chi_table(q); f = dual(q)
    # c(a) = sum_t f(t) chi(a-t) via exact integer FFT-free circular convolution (numpy int64, O(q^2) memory-light loop)
    T = np.nonzero(f)[0]; c = np.zeros(q, dtype=np.int64); a = np.arange(q)
    for t in T: c += f[t]*ch[(a - t) % q]
    N = np.maximum(np.abs(f), np.abs(c)).sum()
    log(f"q={q}: dual pair F={f.sum()} l1={np.abs(f).sum()} sigma(formula)={N-q} ((q-1)/2={(q-1)//2})")
open('check4_out.txt', 'w').write('\n'.join(out) + '\n')
