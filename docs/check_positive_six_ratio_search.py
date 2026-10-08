"""Exact bounded search for the central-height ratio.

This is a finite diagnostic on integer ratio tuples 1 <= t_1 < ... < t_5
<= 30.  The central vector is computed by the closed barycentric formula;
the few tuples near the best ratio are then checked by full Gaussian
normalization and all-edge radius arithmetic.
"""

from fractions import Fraction
from itertools import combinations
from math import gcd, lcm

from check_positive_six_cotangent_search import inspect


def scale_for(ts):
    scale = 1
    for x, y in combinations(ts, 2):
        d = abs(x - y)
        scale = lcm(scale, d // gcd(d, x * y + 1))
    return scale


def central_closed(ts):
    """Primitive central vector, using r_0=1 and r_i=(t_i-i)/(t_i+i)."""
    values = [Fraction(1)]
    for i, t in enumerate(ts):
        denominator = 1
        for j, s in enumerate(ts):
            if i != j:
                denominator *= t - s
        values.append(Fraction(-(t * t + 1) ** 2, denominator))
    common = 1
    for value in values:
        common = lcm(common, value.denominator)
    integers = [value.numerator * (common // value.denominator)
                for value in values]
    content = 0
    for value in integers:
        content = gcd(content, abs(value))
    return tuple(value // content for value in integers)


def main():
    candidates = []
    global_best = None
    for ts in combinations(range(1, 31), 5):
        central = central_closed(ts)
        ratio = Fraction(max(map(abs, central)), ts[0])
        if global_best is None or ratio < global_best[0]:
            global_best = (ratio, ts, central)
        if ratio <= Fraction(841, 2):
            candidates.append((ratio, ts, central))

    candidates.sort()
    verified = []
    for ratio, ts, central in candidates:
        data = inspect(ts)
        actual = tuple(data["central"])
        assert actual == central or actual == tuple(-x for x in central)
        verified.append((ratio, ts, data))

    collision_free = [row for row in verified if not row[2]["collisions"]]
    assert collision_free
    best = min(collision_free)
    assert all(row[2]["collisions"] for row in verified if row[0] < best[0])
    print(f"grid size: {len(list(combinations(range(1, 31), 5)))}")
    print(f"unfiltered minimum ratio: {global_best}")
    print(f"tuples at or below first collision-free ratio: {len(candidates)}")
    print(f"all lower-ratio tuples have collisions: {best[0]}")
    print(best)


if __name__ == "__main__":
    main()
