"""Checks of the ingredients of (V.3) and of h <= Q/(2s)."""
import random, itertools, sys
from f4poly import *

random.seed(31)


def deg_lin(T, v):
    return deg(linfac(T, v))


def main():
    # 1. norm identity and raw difference formula, exact polynomial arithmetic
    for d in range(3, 7):
        Q = 1 << d
        fs = forms_and_sides(d)
        V = range(Q)
        for T in [0, 1, 5 % Q, (1 << d), (1 << (d + 1)) ^ 3, (1 << (d + 3)) ^ 6, random.getrandbits(d + 4)]:
            PV = prod(linfac(T, v) for v in V)
            assert all(x in (0, 1) for x in PV), "norm not in F_2[t]"
            Fs = {}
            for (w, ker, co) in fs:
                for side in (ker, co):
                    F = prod(linfac(T, v) for v in side)
                    assert mul(F, conj(F)) == PV
                    Fs[side] = F
            npair = 0
            pairs = list(itertools.combinations(range(len(fs)), 2))
            if len(pairs) > 60:
                pairs = random.sample(pairs, 60)
            for i, j in pairs:
                wi, ki, ci = fs[i]
                wj, kj, cj = fs[j]
                J = ki & kj
                assert len(J) == Q // 4
                for Si in (ki, ci):
                    for Sj in (kj, cj):
                        a = next(iter(Si & Sj))
                        PJ1 = prod(from_f2(1 ^ v) for v in J)
                        PJY = prod(linfac(T, a ^ v) for v in J)
                        assert add(Fs[Si], Fs[Sj]) == mul(PJ1, PJY), (d, T, i, j)
                        npair += 1
            print(f"d={d} T={T:#x}: norm identity and F_i-F_j=P_J(1)P_J(Y+a) exact on {npair} side pairs")

    # 2. deg P_J(1) >= (Q/4)(d-4), exhaustive over pairs of forms, d <= 9
    for d in range(3, 10):
        Q = 1 << d
        odd = list(range(1, Q, 2))
        kers = {w: [v for v in range(Q) if parity(v & w) == 0] for w in odd}
        mn = None
        for wi, wj in itertools.combinations(odd, 2):
            J = set(kers[wi]) & set(kers[wj])
            s = sum(f2deg(1 ^ v) for v in J)
            mn = s if mn is None else min(mn, s)
        mx_formula = (d - 2) * 2 ** (d - 2)
        print(f"d={d}: min deg P_J(1) = {mn}, claimed lower bound (Q/4)(d-4) = {Q//4*(d-4)}, "
              f"exact sharper bound (Q/4)(d-4)+2 = {Q//4*(d-4)+2}; ok={mn >= Q//4*(d-4)}")
        assert mn >= Q // 4 * (d - 4)

    # 3. sum_v (B - deg(Y+v)) <= 2Q, exhaustive over all T of degree < d+2 (incl. 0)
    for d in range(2, 11):
        Q = 1 << d
        worst = 0
        for T in range(1 << (d + 2)):
            B = max(d, f2deg(T) if T else -1)
            s = sum(B - deg_lin(T, v) for v in range(Q))
            worst = max(worst, s)
        print(f"d={d}: max_T sum_v (B-deg(Y+v)) = {worst}  (bound 2Q = {2*Q}; exact 2^(d+1)-2 = {2*Q-2})")
        assert worst <= 2 * Q

    # 4. h = |cap S_i| <= Q/(2s) on random selections, d = 4..8
    for d in range(4, 9):
        Q = 1 << d
        fs = forms_and_sides(d)
        bad = 0
        for _ in range(3000):
            s = random.randint(2, len(fs))
            sel = random.sample(fs, s)
            C = set(range(Q))
            for (w, k, c) in sel:
                C &= (k if random.random() < 0.5 else c)
            if len(C) * 2 * s > Q:
                bad += 1
        print(f"d={d}: h<=Q/(2s) violated in {bad} of 3000 random selections")
        assert bad == 0


main()
