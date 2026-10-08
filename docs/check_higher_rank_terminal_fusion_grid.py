"""Exact finite SL_n fusion/coherent determinant comparison.

Checks n=2..8, d=2..12 by default. This does not verify the separate
Gaussian residual coefficient-height premise.
"""

from collections import defaultdict
from fractions import Fraction
from math import comb, factorial


def alcove_states(n, level):
    def rec(prefix, remaining):
        if len(prefix) == n - 1:
            yield tuple(prefix)
            return
        for value in range(remaining + 1):
            yield from rec(prefix + [value], remaining - value)
    return tuple(rec([], level))


def fundamental_steps(state, n, level):
    for row in range(n):
        target = list(state)
        if row == 0:
            target[0] += 1
        elif row == n - 1:
            target[n - 2] -= 1
        else:
            target[row - 1] -= 1
            target[row] += 1
        if min(target) >= 0 and sum(target) <= level:
            yield tuple(target)


def casimir_times_n(state, n):
    rows = [sum(state[i:]) for i in range(n - 1)] + [0]
    total = sum(rows)
    return n * sum(rows[i] * (rows[i] + n - 1 - 2 * i)
                   for i in range(n)) - total * total


def source_dimension(n, d):
    numerator = factorial(n * d)
    denominator = 1
    for i in range(n):
        numerator *= factorial(i)
        denominator *= factorial(d + i)
    assert numerator % denominator == 0
    return numerator // denominator


def case(n, d):
    m = n * d
    level = d - 1
    kappa = d + n - 1
    states = alcove_states(n, level)
    transitions = {state: tuple(fundamental_steps(state, n, level))
                   for state in states}
    zero = (0,) * (n - 1)
    history = [{zero: 1}]
    for _ in range(m):
        next_counts = defaultdict(int)
        for state, count in history[-1].items():
            for target in transitions[state]:
                next_counts[target] += count
        history.append(dict(next_counts))
    h = history[m][zero]
    source = source_dimension(n, d)
    assert 0 < h < source

    denominator = 2 * n * kappa
    fundamental_casimir_times_n = n * n - 1

    def tau_numerator(s):
        total = -s * h * fundamental_casimir_times_n
        other = history[m - s]
        for state, multiplicity in history[s].items():
            total += (casimir_times_n(state, n) * multiplicity
                      * other.get(tuple(reversed(state)), 0))
        return total

    tau_pair_num = tau_numerator(2)
    tau_pair = Fraction(tau_pair_num, denominator)
    delta = (-Fraction(h * fundamental_casimir_times_n, n * kappa)
             - (m - 1) * tau_pair)
    assert delta.denominator == 1 and delta > 0
    delta = int(delta)

    nu = [0] * (m + 1)
    for s in range(2, m // 2 + 1):
        order = Fraction(tau_numerator(s)
                         - comb(s, 2) * tau_pair_num, denominator)
        assert order.denominator == 1 and order >= 0
        nu[s] = int(order)
    for s in range(m // 2 + 1, m + 1):
        order = Fraction(delta * (2 * s - m), 2) + nu[m - s]
        assert order.denominator == 1 and order >= 0
        nu[s] = int(order)

    B = m * delta * (1 << (m - 3)) - sum(
        comb(m, s) * nu[s] for s in range(m + 1))
    single_pair_threshold = 1 << ((n - 2) * d + 1)
    threshold = (1 << m) - 2 * ((1 << d) - 1) ** n + ((1 << d) - 2) ** n
    assert threshold >= single_pair_threshold
    assert B >= h * threshold
    return dict(n=n, d=d, m=m, source=source, h=h, delta=delta,
                B=B, threshold=threshold, nu=tuple(nu))


def main():
    results = {}
    for n in range(2, 9):
        for d in range(2, 13):
            results[n, d] = case(n, d)
    equalities = [(n, d) for (n, d), item in results.items()
                  if item["B"] == item["h"] * item["threshold"]]
    assert equalities == [(2, 2), (3, 2)]
    for n in range(2, 9):
        for d in range(2, 9):
            first = results[n, d]
            second = results[d, n]
            assert first["source"] == second["source"]
            assert first["h"] + second["h"] == first["source"]
            assert first["delta"] == second["delta"]
            assert first["B"] == second["B"]
    for n, d in ((2, 2), (3, 2), (4, 2), (3, 4), (4, 3), (6, 2)):
        item = results[n, d]
        ratio = item["B"] / (item["h"] * item["threshold"])
        print(f"(n,d)=({n},{d}): h={item['h']}, delta={item['delta']}, "
              f"B={item['B']}, ratio={ratio:.8f}")
    print("77 exact cases: all-pair floor met only at (2,2),(3,2); no strict gain.")


if __name__ == "__main__":
    main()
