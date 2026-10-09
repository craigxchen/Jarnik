"""C_bal for S_M-symmetric cut polytopes: compare the rigorous peeling upper bound for reg(S(nU))
with (2 - 1/ceil(M/2)) * w_bal (w_bal = exact width of balanced cuts) and with 2*adv."""
import sys, random, math
from fractions import Fraction as Fr
from peeldp import bound
from clib import sym_adversary
random.seed(int(sys.argv[4]))
M = int(sys.argv[1]); trials = int(sys.argv[2]); R = int(sys.argv[3]); n = int(sys.argv[5])
k = math.ceil(M / 2); q = M // 2
def tops(v):
    vs = sorted(v, reverse=True); out = []; a = 0
    for b in range(1, M):
        a += vs[b-1]; out.append(a)
    return out
worst = None; worst2 = None
for t in range(trials):
    s = random.randint(-M + 1, M - 1)
    gens = []
    for _ in range(random.randint(1, 3)):
        v = [random.randint(-R, R) for _ in range(M - 1)]; v.append(s - sum(v)); gens.append(v)
    U = [max(tops(v)[b] for v in gens) for b in range(M - 1)]
    w = {}
    for b in range(1, q + 1):
        w[b] = max(tops(v)[b-1] for v in gens) - min(s - tops(v)[M-b-1] for v in gens)
    if w[q] == 0: continue
    adv, _ = sym_adversary(w, M)
    pb = bound(M, tuple(n * u for u in U), n * s)
    r1 = Fr(pb) / ((2 - Fr(1, k)) * n * w[q])
    r2 = Fr(pb) / (2 * n * adv)
    if worst is None or r1 > worst[0]: worst = (r1, pb, w, U, s)
    if worst2 is None or r2 > worst2[0]: worst2 = (r2, pb, w, U, s)
print(f"M={M} n={n}: max peel/((2-1/k) n w_bal) = {float(worst[0]):.4f} {worst[1:]} ; max peel/(2 n adv) = {float(worst2[0]):.4f}")
