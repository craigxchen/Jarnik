"""Exact checks for the two erasure estimates used in Part VIII.

Balanced vectors of l1 norm two and four in the Paley tests are covered up
to affine relabelling of the field and negating the entire vector. The
Newton test uses all disjoint positive and negative multisets, with their
multiplicities. Python standard library only.
"""
from itertools import combinations_with_replacement
from math import comb, isqrt


def paley_rows(q):
    squares = {x * x % q for x in range(1, q)}
    return [[0 if x == j else (1 if (j - x) % q in squares else -1)
             for x in range(q)] for j in range(q)]


def primitive_shapes(q, length):
    if length == 2:
        yield [0], [1]
        return
    # If one side is doubled, make it the positive side, translate it to
    # zero and dilate to send the first negative point to one.
    for j in range(2, q):
        yield [0, 0], [1, j]
    # Two distinct points on each side. Send the positive points to zero
    # and one in this order. Doubling both sides would not be primitive.
    for j in range(2, q):
        for k in range(j + 1, q):
            yield [0, 1], [j, k]


def check_paley(q):
    rows = paley_rows(q)
    for h in (2, 4):
        minimum = {n: q for n in (2, 3, 4, 5, 7)}
        for positive, negative in primitive_shapes(q, h):
            w = [sum(rows[j][x] for j in positive)
                 - sum(rows[j][x] for j in negative) for x in range(q)]
            for n in minimum:
                minimum[n] = min(minimum[n], sum(a % n != 0 for a in w))
        threshold = isqrt(q - 1) + 1 - h
        for n in (3, 4, 5, 7):
            assert minimum[n] >= threshold
        if h == 2:
            assert minimum[2] == 2 < threshold  # mod 2 cannot replace mod 4
        print(f"q={q}, length={h}, required={threshold}, minima={minimum}")


def check_newton(q, t):
    multisets = list(combinations_with_replacement(range(q), t))
    supports = [set(a) for a in multisets]
    power_sums = [tuple(sum(pow(x, s, q) for x in a) % q
                        for s in range(1, t + 1)) for a in multisets]
    eligible = [0] * (t + 1)
    for i, a in enumerate(power_sums):
        for j, b in enumerate(power_sums):
            if not supports[i].isdisjoint(supports[j]):
                continue
            eligible[0] += 1
            for k, (x, y) in enumerate(zip(a, b), 1):
                if x != y:
                    break
                eligible[k] += 1
    assert eligible[-1] == 0
    for k, count in enumerate(eligible):
        bound = comb(q + t - 1, t) * comb(q + t - k - 1, t - k)
        assert count <= bound
    print(f"Newton q={q}, size of each multiset={t}, "
          f"counts matching k=0,...,{t} moments={eligible}")


if __name__ == "__main__":
    for q in (43, 127):
        check_paley(q)
    for t in range(1, 5):
        check_newton(7, t)
    check_newton(11, 3)
