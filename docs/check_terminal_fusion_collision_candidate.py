"""Exact checks for terminal fusion collision content and its height table.

The source bridge and pair subtraction are proved in
terminal_kernel_fusion_content.md. This checker verifies the finite fusion
recurrence and arithmetic consequences; it does not prove the KZ-connection
or determinant-frame arguments in that note.
"""

from fractions import Fraction
from math import comb, factorial


def fusion_histories(d):
    level = d - 1
    histories = [{(0, 0): 1}]
    for _ in range(3 * d):
        nxt = {}
        for (a, b), multiplicity in histories[-1].items():
            for aa, bb in ((a + 1, b), (a - 1, b + 1), (a, b - 1)):
                if aa >= 0 and bb >= 0 and aa + bb <= level:
                    nxt[(aa, bb)] = nxt.get((aa, bb), 0) + multiplicity
        histories.append(nxt)
    return histories


def joint_kernel_dimension(n, degree):
    """Bounded-walk formula for h(n, degree), with L=n-2*degree."""
    width = n - 2 * degree
    if width < 0:
        return 0
    answer = 0
    period = width + 2
    for j in range(-n - 1, n + 2):
        top = degree + j * period
        bottom = degree - 1 + j * period
        if 0 <= top <= n:
            answer += comb(n, top)
        if 0 <= bottom <= n:
            answer -= comb(n, bottom)
    return answer


def casimir(a, b):
    return Fraction(2 * (a * a + a * b + b * b), 3) + 2 * (a + b)


def candidate(k):
    d = 2 * k
    m = 3 * d
    histories = fusion_histories(d)
    h = histories[m].get((0, 0), 0)
    source_dimension = (2 * factorial(m)
                        // (factorial(d) * factorial(d + 1)
                            * factorial(d + 2)))
    coherent_rank = joint_kernel_dimension(m, d)
    assert h == source_dimension - coherent_rank
    tau = []
    for r in range(m + 1):
        weighted = sum(
            (casimir(a, b) * mult * histories[m - r].get((b, a), 0)
             for (a, b), mult in histories[r].items()),
            Fraction(0),
        )
        tau.append((weighted - Fraction(r * h * 8, 3))
                    / (2 * (d + 2)))

    # The pair subtraction is the numerator order proved in
    # terminal_kernel_fusion_content.md, Section 3.
    tau2 = tau[2]
    nu = [tau[r] - comb(r, 2) * tau2 for r in range(m + 1)]
    delta = joint_kernel_dimension(m - 1, d)
    excess = [nu[r] - delta * max(0, r - m // 2)
              for r in range(m + 1)]
    total_content = sum(comb(m, r) * nu[r] for r in range(m + 1))
    covolume_numerator = m * delta * 2 ** (m - 3) - total_content
    assert all(value.denominator == 1 and value >= 0 for value in nu)
    assert all(nu[r] - nu[m - r] == delta * (r - m // 2)
               for r in range(m + 1))
    if k == 1:
        assert nu == [0, 0, 0, 0, 1, 2, 3]
    return {
        "k": k,
        "d": d,
        "m": m,
        "h": h,
        "delta": delta,
        "tau2": tau2,
        "nu": nu,
        "excess": excess,
        "total_content": total_content,
        "covolume_numerator": covolume_numerator,
        "covolume_over_h": Fraction(covolume_numerator, h),
        "threshold": 2 ** (d + 1),
    }


def main():
    print("k | h | delta | B/h | 2^(d+1) | nonzero excess positions")
    for k in range(1, 16):
        result = candidate(k)
        assert result["covolume_over_h"] > result["threshold"]
        nonzero_excess = [(r, value) for r, value in enumerate(result["excess"])
                          if value]
        print(
            f'{k} | {result["h"]} | {result["delta"]} | '
            f'{result["covolume_over_h"]} | {result["threshold"]} | '
            f'{nonzero_excess}'
        )
    print("PASS: fusion dimensions, collision orders, and height data through k=15.")


if __name__ == "__main__":
    main()
