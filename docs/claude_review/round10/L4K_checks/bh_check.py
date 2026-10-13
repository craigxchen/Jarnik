# Exact check of Theorem 1(b): a_j = 2(Lam + (K+1) beta_j) + nu_j, Lam = K(K+1)B + K, beta a B_h set (h = K//2)
# in [1,B]: for every parity pattern nu, every permutation-free choice, every v with 1<=||v||_1<=K,
# <v,a> != 0 as an integer and |<v,a>| <= A_K = K(2(K+1)^2 B + 2K + 1).
import itertools
from lib import monos

def is_Bh(beta, h):
    seen = set()
    for comb in itertools.combinations_with_replacement(beta, h):
        s = sum(comb)
        if s in seen: return False
        seen.add(s)
    return True

def mian_chowla_Bh(r, h):
    beta = []
    x = 1
    while len(beta) < r:
        if is_Bh(beta + [x], h): beta.append(x)
        x += 1
    return beta

allok = True
for K in range(1, 6):
    h = K // 2
    for r in range(2, 6 if K >= 4 else 7):
        beta = list(range(1, r + 1)) if h <= 1 else mian_chowla_Bh(r, h)
        assert h == 0 or is_Bh(beta, h)
        B = max(beta); Lam = K * (K + 1) * B + K
        AK = K * (2 * (K + 1) ** 2 * B + 2 * K + 1)
        V = monos(r, K)
        worst = 0; ok = True
        for nuv in itertools.product((0, 1), repeat=r):
            for perm in ([tuple(range(r))] if r > 4 else itertools.permutations(range(r))):
                a = [2 * (Lam + (K + 1) * beta[perm[j]]) + nuv[j] for j in range(r)]
                for v in V:
                    s = sum(x * y for x, y in zip(v, a))
                    if s == 0: ok = False
                    worst = max(worst, abs(s))
        assert worst <= AK
        allok &= ok
        print(f"K={K} h={h} r={r} beta={beta} B={B} A_K={AK} max|<v,a>|={worst} nonzero_all={ok}")
print("ALL OK" if allok else "FAIL")
