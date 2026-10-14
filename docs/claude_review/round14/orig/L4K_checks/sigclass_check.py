# Of the label vectors with the genuine collision signature at q, how many are in the genuine class?
import math, numpy as np, sys
from forced_sig import sig_table
from forced_coh import load_set
K, r, tag = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
P, F = load_set(K, r, tag)
for q in F:
    t = sig_table(P, K, q); n = t['n']; Lc = t['L'] % n
    units = [u for u in range(1, n + 1) if math.gcd(u, n) == 1] if n > 1 else [1]
    G = set(tuple(int(x) for x in (u * Lc) % n) for u in units)
    idx = np.where(t['inv'] == t['genuine_key'])[0]
    ing = sum(tuple(int(x) for x in t['ells'][i] % n) in G for i in idx)
    print(f"K={K} r={r} {tag} q={q} n={n}: vectors with genuine signature {len(idx)}, in genuine class {ing}")
