"""Referee check of sharp.md Prop. 4.3 (= lower.md Prop. 3.3, canonical assignment) and Lemma 4.1.

Own construction of the canonical one-flip Paley profile with b copies:
  labels a in F_q, physical columns (a, k), k < b; row x (finite) flipped in (x, 0) with sibling
  (x, 1); row infinity flipped in (a*, 2) with sibling (a*, 3).
Checks (exact integers):
  1. twins: S_f(x) - S_s(x) = -2 H(x, alpha(x)) e_x, rank S = M, every column nonconstant;
  2. identity (3.1): 2B = b F + sum_x (|T_alpha(x) - 2 c_x H(x,alpha(x))| - |T_alpha(x)|);
  3. 2B - r >= b - 4 for all pairs (exhaustive), >= b + 2M - 8 for all signed columns,
     >= b + M - 16 for all signed half-sums, >= (b/16 - 2) M + b for other ternary c,
     >= h M (b-4)/2 + b for amplitude h >= 2;  the "kappa_b M" summary line;
  4. slack inequality s_c >= (tau/2)(2B - r) - 1/2 with random weights w_j in [tau, tau + 1/r];
  5. adversarial search for the smallest 2B - r over non-pair characters (hill climbing on 2B).
"""
import itertools, math, random
import numpy as np

rng = random.Random(77)


def chi(a, p):
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def paley(q):
    M = q + 1
    H = np.empty((M, q), dtype=np.int64)
    for x in range(q):
        for a in range(q):
            H[x, a] = -chi(a - x, q) - (1 if a == x else 0)
    H[q, :] = 1
    return H


def profile(q, b, astar=5):
    H = paley(q)
    M = q + 1
    r = b * q
    col = lambda a, k: a * b + k
    S = np.repeat(H, b, axis=1)
    f, s, al = {}, {}, {}
    for x in range(q):
        f[x], s[x], al[x] = col(x, 0), col(x, 1), x
    f[q], s[q], al[q] = col(astar, 2), col(astar, 3), astar
    for x in range(M):
        S[x, f[x]] *= -1
    return H, S, f, s, al


ok = True
lines = []
for q in (43, 47, 59, 67):
    for b in (33, 40):
        H, S, f, s, al = profile(q, b)
        M, r = S.shape
        # 1. twins / rank / nonconstant
        tw = all(((S[:, f[x]] - S[:, s[x]]) == -2 * H[x, al[x]] * np.eye(M, dtype=np.int64)[x]).all() for x in range(M))
        rk = np.linalg.matrix_rank(S.astype(float)) == M
        nc = all(len(set(S[:, j])) == 2 for j in range(r))
        ok &= tw and rk and nc

        def twoB(c):
            return int(np.abs(c @ S).sum())

        def ident(c):
            T = c @ H
            F = int(np.abs(T).sum())
            fl = sum(abs(int(T[al[x]]) - 2 * int(c[x]) * int(H[x, al[x]])) - abs(int(T[al[x]])) for x in range(M))
            return b * F + fl

        # pairs exhaustive
        mp = 10 ** 9
        idok = True
        for x in range(M):
            for y in range(x + 1, M):
                c = np.zeros(M, dtype=np.int64)
                c[x], c[y] = 1, -1
                v = twoB(c)
                idok &= v == ident(c)
                mp = min(mp, v - r)
        ok &= mp >= b - 4 and idok
        # columns
        mc = min(twoB(sg * H[:, a]) - r for a in range(q) for sg in (1, -1))
        ok &= mc >= b + 2 * M - 8
        # half-sums
        mh = 10 ** 9
        for a, a2 in itertools.combinations(range(q), 2):
            for s1, s2 in ((1, 1), (1, -1), (-1, 1), (-1, -1)):
                g = (s1 * H[:, a] + s2 * H[:, a2]) // 2
                mh = min(mh, twoB(g) - r)
        ok &= mh >= b + M - 16
        # other ternary random + amplitude
        mo = 10 ** 9
        ma_ratio = 10 ** 9
        for it in range(3000):
            h = 2 * rng.randrange(2, M // 2 + 1)
            Sx = rng.sample(range(M), h)
            c = np.zeros(M, dtype=np.int64)
            for i, x in enumerate(Sx):
                c[x] = 1 if i < h // 2 else -1
            v = twoB(c)
            idok &= v == ident(c)
            Tabs = np.abs(c @ H)
            exc = ((Tabs == M).sum() == 1 and (Tabs == 0).sum() == q - 1) or ((Tabs == M // 2).sum() == 2 and (Tabs == 0).sum() == q - 2)
            if not exc:
                mo = min(mo, v - r)
        for it in range(1500):
            hh = rng.choice([2, 3, 5])
            c = np.array([rng.randint(-hh, hh) for _ in range(M)], dtype=np.int64)
            c[rng.randrange(M)] -= c.sum()
            h = int(np.abs(c).max())
            if h < 2:
                continue
            v = twoB(c)
            idok &= v == ident(c)
            bound = h * M * (b - 4) / 2 + b
            ma_ratio = min(ma_ratio, (v - r) / bound)
        ok &= mo >= (b / 16 - 2) * M + b and ma_ratio >= 1 and idok
        kb = min(1, b / 16 - 2)
        allother = min(mc, mh, mo)
        ok &= allother >= kb * M
        # 4. slack inequality with random weights
        tau = 50.0
        worst_sl = 10 ** 9
        for it in range(300):
            w = tau + rng.random() * (1.0 / r) * np.ones(r)
            w = tau + np.array([rng.random() / r for _ in range(r)])
            h = 2 * rng.randrange(1, M // 2 + 1)
            Sx = rng.sample(range(M), h)
            c = np.zeros(M, dtype=np.int64)
            for i, x in enumerate(Sx):
                c[x] = 1 if i < h // 2 else -1
            tauj = c @ S
            sc = 0.5 * float((w * (np.abs(tauj) - 1)).sum())
            lb = (tau / 2) * (twoB(c) - r) - 0.5
            worst_sl = min(worst_sl, sc - lb)
        ok &= worst_sl >= -1e-9
        # 5. adversarial hill climb for small 2B - r among non-pairs (ternary)
        best = 10 ** 9
        for start in range(25):
            h = 2 * rng.randrange(2, M // 2 + 1)
            Sx = rng.sample(range(M), h)
            c = np.zeros(M, dtype=np.int64)
            for i, x in enumerate(Sx):
                c[x] = 1 if i < h // 2 else -1
            cur = twoB(c)
            for step in range(1500):
                c2 = c.copy()
                x, y = rng.sample(range(M), 2)
                c2[x] += 1
                c2[y] -= 1
                if np.abs(c2).max() > 1 or np.abs(c2).sum() <= 2:
                    continue
                v2 = twoB(c2)
                if v2 <= cur:
                    c, cur = c2, v2
            best = min(best, cur - r)
        ok &= best >= kb * M
        lines.append(f"q={q:3d} b={b}: twins/rank/nonconst {tw and rk and nc}; identity(3.1) {idok}; min 2B-r: pairs {mp} (b-4={b-4}), "
                     f"columns {mc} (b+2M-8={b+2*M-8}), half-sums {mh} (b+M-16={b+M-16}), other ternary (random) {mo} "
                     f"(proved {(b/16-2)*M+b:.1f}), hill-climb non-pair {best}; amplitude min ratio to bound {ma_ratio:.2f}; "
                     f"kappa_b M = {kb*M:.2f}; slack ineq min gap {worst_sl:.3f}")
lines.append("ALL REFEREE SLACK CHECKS PASSED" if ok else "FAILURE")
print("\n".join(lines))
