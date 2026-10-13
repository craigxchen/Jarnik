"""Independent check of Theorem 1(b): a_j = 2(Lambda+(K+1)beta_sigma(j)) + nu_j has <v,a> != 0 and
|<v,a>| <= A_K for all v in V_K, every parity pattern nu, random permutations sigma."""
import itertools, random
import numpy as np
from rlib import VK

def greedy_Bh(r, h):
    if h <= 1: return list(range(1, r + 1))
    S = []
    x = 0
    while len(S) < r:
        x += 1
        T = S + [x]
        sums = set(); ok = True
        for c in itertools.combinations_with_replacement(T, h):
            s = sum(c)
            if s in sums: ok = False; break
            sums.add(s)
        if ok: S = T
    return S

rng = random.Random(7)
allok = True
for K in range(1, 7):
    h = K // 2
    for r in range(2, 7 if K <= 4 else 5):
        beta = greedy_Bh(r, h); B = max(beta)
        Lam = K * (K + 1) * B + K
        AK = K * (2 * (K + 1) ** 2 * B + 2 * K + 1)
        V = np.array(VK(r, K))
        mx = 0; zero = 0
        pats = list(itertools.product([0, 1], repeat=r))
        for nu in pats:
            for _ in range(6):
                sig = list(range(r)); rng.shuffle(sig)
                a = np.array([2 * (Lam + (K + 1) * beta[sig[j]]) + nu[j] for j in range(r)])
                t = V @ a
                zero += int((t == 0).sum()); mx = max(mx, int(np.abs(t).max()))
        ok = (zero == 0 and mx <= AK)
        allok &= ok
        print(f"K={K} r={r} h={h} B={B} A_K={AK} max|<v,a>|={mx} zeros={zero} {'OK' if ok else 'FAIL'}")
print("ALL OK" if allok else "FAILURE")
