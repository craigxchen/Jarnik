#!/usr/bin/env python3
"""R1. Independent evaluation of F*(s) (two_thirds Thm 5.2) and of the inputs (T2) used by stability.md.

(a) psi8_dp == psi8_brute (all compositions) for every split p and 2 <= s <= 15;
(b) 8F*(s) from rlib agrees with stab_lib.Fstar (the write-up's evaluator) for 2 <= s <= 80 (to 1e-40);
(c) (T2): 8F*(s)/s > 3 s log s - 13.9092 s + 3 + 2 log 2 for 2 <= s <= SMAX (stronger form, with the
    +3+2log2 the write-up quotes), and F*(s) >= binom(s,2) log2 / 2;
(d) (3.0): for M >= 4 and ceil(M/2) <= s <= M: 8F*(s)/s^2 > 3 log(M/2) - K0  (i.e. the second form).
Usage: python3 r1_fstar.py [SMAX]
"""
import sys
import mpmath
sys.path.insert(0, '../stability_checks')
import stab_lib                                   # noqa: E402  (only for the cross-check (b))
from rlib import psi8_dp, psi8_brute, F8vec, sieve  # noqa: E402

mpmath.mp.dps = 50


def F8mp(s):
    return sum(n * mpmath.log(q) for q, n in F8vec(s).items())


def main():
    SMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 250
    ok = True
    # (a)
    cnt = 0
    for s in range(2, 16):
        for p in sieve(s + 2):
            if p % 4 == 1 and p - 1 < s:
                a, b = psi8_dp(s, p), psi8_brute(s, p)
                cnt += 1
                if a != b:
                    ok = False
                    print("psi mismatch", s, p, a, b)
    print(f"(a) DP == brute force for Psi_p on {cnt} (s, p) pairs, s <= 15")
    # (b)
    worst = 0
    for s in range(2, 81):
        d = abs(F8mp(s) / 8 - stab_lib.Fstar(s))
        worst = max(worst, d)
    print(f"(b) max |F*_referee - F*_stab_lib| over 2 <= s <= 80: {mpmath.nstr(worst, 5)}")
    if worst > mpmath.mpf(10) ** -40:
        ok = False
    # (c), (d)
    K0 = mpmath.mpf('13.9092')
    minslack = None
    for s in range(2, SMAX + 1):
        f8 = F8mp(s)
        lhs = f8 / s
        rhs = 3 * s * mpmath.log(s) - K0 * s + 3 + 2 * mpmath.log(2)
        if not lhs > rhs:
            ok = False
            print("T2 fails at", s)
        sl = (lhs - rhs) / s
        if minslack is None or sl < minslack[0]:
            minslack = (sl, s)
        if not f8 / 8 >= mpmath.mpf(s * (s - 1)) / 4 * mpmath.log(2):
            ok = False
    print(f"(c) (T2) with +3+2log2 holds for 2 <= s <= {SMAX}; min (lhs-rhs)/s = "
          f"{mpmath.nstr(minslack[0], 6)} at s = {minslack[1]}")
    worst = None
    for M in range(4, SMAX + 1):
        for s in range((M + 1) // 2, M + 1):
            v = F8mp(s) / s ** 2 - (3 * mpmath.log(mpmath.mpf(M) / 2) - K0)
            if v <= 0:
                ok = False
            if worst is None or v < worst[0]:
                worst = (v, M, s)
    print(f"(d) 8F*(s)/s^2 > 3 log(M/2) - K0 for all ceil(M/2) <= s <= M <= {SMAX}; min margin "
          f"{mpmath.nstr(worst[0], 6)} at (M, s) = {worst[1:]}")
    # tabulate a few values for the report
    for s in (2, 3, 4, 6, 8, 12, 16, 24, 48):
        print(f"   s = {s:3d}: F*(s) = {mpmath.nstr(F8mp(s) / 8, 10)},  8F*/s^2 = {mpmath.nstr(F8mp(s) / s**2, 8)},"
              f"  3 log s - 8F*/s^2 = {mpmath.nstr(3 * mpmath.log(s) - F8mp(s) / s**2, 6)}")
    print("ALL CHECKS PASSED:", ok)


if __name__ == "__main__":
    main()
