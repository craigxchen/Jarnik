"""Exact bounded scan of lattice circles for repeated pair squareclasses.

For every odd squared radius N <= LIMIT, all integer points x^2+y^2=N are
enumerated.  Pair norms use the exact identity d=N/gcd(Re(z*conj(w)),
Im(z*conj(w))).  Only the angular ordering and the displayed C value use
ordinary floating point.
"""

from itertools import combinations
from math import atan2, gcd, isqrt, pi


LIMIT = 100_000


def points(N):
    out = []
    for x in range(isqrt(N) + 1):
        y = isqrt(N - x * x)
        if y * y != N - x * x:
            continue
        xs = (x, -x) if x else (0,)
        ys = (y, -y) if y else (0,)
        out.extend((a, b) for a in xs for b in ys)
    return list(dict.fromkeys(out))


def squareclass(value):
    out = []
    p = 2
    while p * p <= value:
        exponent = 0
        while value % p == 0:
            value //= p
            exponent += 1
        if exponent & 1:
            out.append(p)
        p += 1
    if value > 1:
        out.append(value)
    return tuple(out)


def pair_data(N, left, right):
    real = left[0] * right[0] + left[1] * right[1]
    imag = left[1] * right[0] - left[0] * right[1]
    content = gcd(abs(real), abs(imag))
    assert content and N % content == 0
    d = N // content
    return d, squareclass(d)


def scan_circle(N):
    ordered = sorted(points(N), key=lambda z: atan2(z[1], z[0]) % (2 * pi))
    if len(ordered) < 4:
        return None
    doubled = ordered + ordered
    best = None
    for start in range(len(ordered)):
        block = doubled[start:start + 4]
        first = atan2(block[0][1], block[0][0]) % (2 * pi)
        last = atan2(block[-1][1], block[-1][0]) % (2 * pi)
        width = last - first
        if start + 4 > len(ordered):
            width += 2 * pi
        C = N ** 0.25 * width
        edges = []
        for i, j in combinations(range(4), 2):
            d, label = pair_data(N, block[i], block[j])
            edges.append((i, j, d, label))
        repeated = [pair for pair in combinations(edges, 2)
                    if pair[0][3] == pair[1][3]
                    and pair[0][2] != pair[1][2]
                    and not set(pair[0][:2]) & set(pair[1][:2])]
        if repeated and (best is None or C < best[0]):
            best = (C, block, edges)
    return best


def main():
    checked = 0
    best = None
    under_two = 0
    for N in range(1, LIMIT + 1, 2):
        result = scan_circle(N)
        if result is None:
            continue
        checked += 1
        if best is None or result[0] < best[0]:
            best = (result[0], N, result[1], result[2])
        if result[0] <= 2:
            under_two += 1
    print(f"PASS: checked {checked} odd circles with at least four lattice points.")
    print(f"PASS: {under_two} had a repeated-label four-window with C<=2.")
    if best is not None:
        C, N, block, edges = best
        print(f"BEST: N={N}, C={C:.12g}, points={block}")
        print(f"      pair data={edges}")


if __name__ == "__main__":
    main()
