# ext_test.py W.pkl K2 : evaluate a GB witness (labels on F_K) against the norm-<=K2 budgets m_(v,s) <= |D_s(Z_v)|.
# Violations caused by the stored labels alone are a LOWER bound for every extension (more labels only add demand).
import sys, pickle, math
from collections import Counter
import slverify as sv
W = pickle.load(open(sys.argv[1], "rb")); K2 = int(sys.argv[2]); P = W["P"]; r = len(P); pis = [sv.gp(p) for p in P]
V = sv.monos(r, K2); mv = {}
for q in W["F"]:
    lab = W["labels"][q]; top = lab["top"]; A = lab["A"]; m = q - (1 if q % 4 == 1 else -1); n = m // 4
    tab = sv.dlog_table(tuple(lab["g1"]), q, m); iota = tab[(0, 1)]; Mtop = m * q ** (A - 1)
    for v in V:
        d = sum(a * b for a, b in zip(v, top)) % Mtop
        if d % n: continue
        s = [u for u in range(4) if (u * iota) % m == d % m][0]; b = 1
        while b < A and d % (n * q ** b) == 0: b += 1
        persist = d % (n * q ** (A - 1)) == 0
        mv[(v, s)] = mv.get((v, s), 1) * q ** b
bad = Counter(); worst = []
for v in V:
    X, Y = sv.Zv(pis, v); D = (Y, X - Y, X, X + Y)
    for s in range(4):
        mm = mv.get((v, s), 1)
        if mm > abs(D[s]): bad[sum(map(abs, v))] += 1; worst.append((math.log(mm / abs(D[s])), v, s, mm, D[s]))
worst.sort(reverse=True)
print(f"{sys.argv[1]} at K={K2}: violated (v,s) by norm (both signs counted): {dict(sorted(bad.items()))}; total {sum(bad.values())} of {4*len(V)}")
for t in worst[:5]: print("   ", t)
