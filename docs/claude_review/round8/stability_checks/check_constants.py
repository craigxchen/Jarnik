#!/usr/bin/env python3
"""Numerical constants of stability.md Sections 3-5 (exact integer E-sums, 60-digit logs).

* (4.3): 2 sum_{l = 3 mod 4} E(M, l+1) log l > (M^2/2) log M - 3.9626 M^2 for 2 <= M <= MMAX, and the
  constant (beta_0 + xi)/2 + a_3 + log 4 < 3.9626 from the two_thirds bounds.
* (T2): 8F*(s)/s > 3 s log s - 13.9092 s for 2 <= s <= 60 (F* exact), as a sanity check of the quoted input.
* Proposition 3.4: (M/2)(1 - M^(-1/3))^2 > M/4 exactly when M >= 40.
"""
import sys
import mpmath
from stab_lib import E, primes_upto, Fstar

mpmath.mp.dps = 60


def main():
    MMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
    sigma_hi = mpmath.mpf('0.7553683')
    beta0 = 1 + sigma_hi + mpmath.log(2) / 2
    xi = mpmath.mpf('2.70189')
    a3 = mpmath.mpf('0.1743774')
    const = (beta0 + xi) / 2 + a3 + mpmath.log(4)
    print("beta_0 <=", mpmath.nstr(beta0, 8), "  (beta_0 + xi)/2 + a_3 + log 4 =", mpmath.nstr(const, 8),
          "< 3.9626:", const < mpmath.mpf('3.9626'))
    ls = [l for l in primes_upto(MMAX + 2) if l % 4 == 3]
    worst = None
    ok = True
    for M in range(2, MMAX + 1):
        s = sum(2 * E(M, l + 1) * mpmath.log(l) for l in ls if l + 1 < M)
        rhs = mpmath.mpf(M) ** 2 / 2 * mpmath.log(M) - mpmath.mpf('3.9626') * M * M
        sl = (s - rhs) / M ** 2
        if s <= rhs:
            ok = False
        if worst is None or sl < worst[0]:
            worst = (float(sl), M)
    print(f"(4.3) holds for 2 <= M <= {MMAX}: {ok};  min (lhs - rhs)/M^2 = {worst[0]:.4f} at M = {worst[1]}")
    okT2 = True
    for s in range(2, 61):
        if not 8 * Fstar(s) / s > 3 * s * mpmath.log(s) - mpmath.mpf('13.9092') * s:
            okT2 = False
    print("(T2) 8F*(s)/s > 3 s log s - 13.9092 s for 2 <= s <= 60:", okT2)
    first = None
    for M in range(8, 1000):
        if M / 2 * (1 - M ** (-1 / 3)) ** 2 > M / 4:
            first = M
            break
    print("Prop 3.4 bound exceeds M/4 from M =", first)
    print("ALL CHECKS PASSED:", ok and okT2 and const < mpmath.mpf('3.9626'))


if __name__ == "__main__":
    main()
