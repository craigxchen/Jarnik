# C2: independent exact check of Proposition D in the construction regime (p_j > X), as in L4 §7.
import random, sys
from glib import *
from consD import *

def run(X, r, seed):
    rng = random.Random(seed)
    ps = set()
    while len(ps) < r:
        p = rng.randrange(10**4, 10**6)
        if p % 4 == 1 and is_prime(p): ps.add(p)
    ps = sorted(ps)
    MOD = 4 * 17 * 31 * 47
    pp, used = [], set()
    for p in ps:
        c = p + MOD
        while not (is_prime(c) and c not in used): c += MOD
        used.add(c); pp.append(c)
    zetas = [two_squares(c) for c in pp]
    S = System(ps, X)
    ell, case, L = S.construct_D(zetas)
    cnt = {1: 0, 2: 0, 3: 0, 4: 0}
    for q in S.qs: cnt[case[q]] += 1
    # parity
    for q in S.qs:
        for j in range(r):
            assert ell[q][j] % 2 == S.nu[q][j], ("parity", q, j)
    res = {}
    for k in (1, 2, 3):
        mons = monomials(r, k)
        fails = 0; fails13 = 0
        for v in mons:
            prof = collision_profile(S, ell, v)
            mv = 1
            for s in range(4): mv *= mod_from(prof[s])
            W = Zv(zetas, v)
            I16 = gpow(W, 16)[1]
            assert I16 != 0
            if I16 % mv != 0: fails += 1
            mv13 = 1
            for s in range(4): mv13 *= mod_from([(q, a) for q, a in prof[s] if case[q] != 4])
            if I16 % mv13 != 0: fails13 += 1
        res[k] = (len(mons), fails, fails13)
    print(f"X={X} r={r} seed={seed} cases={cnt} parity OK; per norm k: (#monomials, #fail m_v|Im(Z^16), #fail on cases1-3 only) = {res}")
    return res

allok = True
for X, r, seed in [(400, 12, 1), (1000, 12, 2), (2000, 12, 3), (3000, 12, 4)]:
    res = run(X, r, seed)
    if res[1][1] or res[2][1] or res[3][2]: allok = False
print("Prop D norm<=2 and cases1-3 norm-3:", "ALL OK" if allok else "FAILURE")
