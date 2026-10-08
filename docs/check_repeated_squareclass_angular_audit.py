"""Bounded exact audit for repeated pair-norm squareclasses.

The Gaussian arithmetic (points, pair norms, and squareclasses) is exact.
Angles are used only to order the exact points and to report the normalized
arc width; the search is a bounded numerical audit, not a proof by sampling.
"""

from itertools import combinations, product
from math import atan2, isqrt, pi, sqrt


def gmul(x, y):
    return x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0]


def gpow(x, n):
    out = (1, 0)
    for _ in range(n):
        out = gmul(out, x)
    return out


def split_prime(p):
    for a in range(1, isqrt(p) + 1):
        b = isqrt(p - a * a)
        if a * a + b * b == p:
            return a, b
    raise ValueError(f"{p} is not split")


def is_prime(p):
    return p >= 2 and all(p % q for q in range(2, isqrt(p) + 1))


def systems_points(factors):
    """Return the exact norm N and all exponent-row points for a system."""
    options = []
    N = 1
    for p, width in factors:
        N *= p ** width
        pi_p = split_prime(p)
        options.append([
            gmul(gpow(pi_p, e), gpow((pi_p[0], -pi_p[1]), width - e))
            for e in range(width + 1)
        ])
    points = []
    for row in product(*[range(len(col)) for col in options]):
        z = (1, 0)
        for col, e in zip(options, row):
            z = gmul(z, col[e])
        points.append((z, row))
    assert all(x * x + y * y == N for (x, y), _ in points)
    return N, points


def pair_data(factors, left, right):
    """Return exact pair norm and its rational squareclass label."""
    d = 1
    label = []
    for (p, _), e, f in zip(factors, left, right):
        gap = abs(e - f)
        d *= p ** gap
        if gap & 1:
            label.append(p)
    return d, tuple(label)


def repeated_edges(factors, rows):
    edges = []
    for i, j in combinations(range(len(rows)), 2):
        d, label = pair_data(factors, rows[i], rows[j])
        edges.append((i, j, d, label))
    repeats = []
    for left, right in combinations(edges, 2):
        # Norm injectivity is an already-proved endpoint consequence.  A
        # repeated squareclass is relevant only when the ordinary norms
        # themselves differ.
        if (left[3] == right[3] and left[2] != right[2]
                and not set(left[:2]) & set(right[:2])):
            repeats.append((left, right))
    return edges, repeats


def cyclic_windows(factors, points, size=4):
    """Yield four consecutive angular points, including windows across 0."""
    ordered = sorted(
        points,
        key=lambda item: atan2(item[0][1], item[0][0]) % (2 * pi),
    )
    doubled = ordered + ordered
    # N is the squared radius, hence Delta <= C*N^(-1/4).
    scale = systems_points(factors)[0] ** 0.25
    for start in range(len(ordered)):
        block = doubled[start:start + size]
        first = atan2(block[0][0][1], block[0][0][0]) % (2 * pi)
        last = atan2(block[-1][0][1], block[-1][0][0]) % (2 * pi)
        width = last - first
        if width < 0:
            width += 2 * pi
        # For a wrapped block, the doubled angular lift is the right width.
        if start + size > len(ordered):
            width = (atan2(block[-1][0][1], block[-1][0][0]) % (2 * pi)
                     + 2 * pi - first)
        yield scale * width, [item[1] for item in block]


def audit_system(factors):
    N, points = systems_points(factors)
    best_repeated = float("inf")
    best_all = float("inf")
    repeated_witness = None
    for C, rows in cyclic_windows(factors, points):
        if C < best_all:
            best_all = C
        _, repeats = repeated_edges(factors, rows)
        if repeats and C < best_repeated:
            best_repeated = C
            repeated_witness = (rows, repeats)
    return N, len(points), best_all, best_repeated, repeated_witness


def selected_C(factors, rows):
    N, points = systems_points(factors)
    by_row = {row: z for z, row in points}
    angles = sorted(atan2(by_row[row][1], by_row[row][0]) for row in rows)
    return N ** 0.25 * (angles[-1] - angles[0])


def small_split_primes(limit):
    return [p for p in range(5, limit + 1)
            if is_prime(p) and p % 4 == 1]


def bounded_systems(limit=100):
    primes = small_split_primes(limit)
    # Exhaust all 3-coordinate systems and all widths <=2.  For four
    # coordinates, exhaust the squarefree layer and one doubled coordinate;
    # this is enough to expose the first repeated-label patterns.
    for k in (3, 4):
        for chosen in combinations(primes, k):
            width_patterns = [(1,) * k]
            if k == 3:
                width_patterns += list(product((1, 2), repeat=k))
            else:
                width_patterns += [
                    tuple(2 if i == j else 1 for i in range(k))
                    for j in range(k)
                ]
            for widths in width_patterns:
                yield tuple(zip(chosen, widths))


def main():
    systems = list(bounded_systems())
    counts = {"systems": 0, "repeated": 0, "under_2": 0}
    best = None
    for factors in systems:
        N, m, all_C, repeated_C, witness = audit_system(factors)
        counts["systems"] += 1
        if repeated_C < float("inf"):
            counts["repeated"] += 1
            candidate = (repeated_C, factors, N, m, all_C, witness)
            if best is None or candidate[:3] < best[:3]:
                best = candidate
            if repeated_C <= 2:
                counts["under_2"] += 1
    print(f"PASS: audited {counts['systems']} exact split-prime allocation systems.")
    print(f"PASS: {counts['repeated']} systems contain disjoint repeated labels;")
    print(f"      {counts['under_2']} have a four-point window with C<=2.")
    if best is not None:
        C, factors, N, m, all_C, witness = best
        print(f"BEST repeated-label window: C={C:.12g}, factors={factors}, N={N}, M={m}.")
        print(f"BEST unrestricted four-window in that system: C={all_C:.12g}.")
        rows, repeats = witness
        for left, right in repeats:
            print(f"  repeated label {left[3]}: d={left[2]} on {left[:2]},"
                  f" d={right[2]} on {right[:2]}; rows={rows}")
    # Include the known exact four-point Pell fixture as an explicit angular
    # regression.  It has C<2 and no repeated pair squareclass.
    pell = ((5, 1), (233, 1), (89, 1), (3461, 1))
    pell_rows = ((1, 1, 1, 1), (0, 1, 0, 0),
                 (0, 0, 1, 0), (0, 0, 0, 1))
    N, _ = systems_points(pell)
    assert selected_C(pell, pell_rows) < 2
    _, repeats = repeated_edges(pell, list(pell_rows))
    assert not repeats
    pell_edges, _ = repeated_edges(pell, list(pell_rows))
    assert len({edge[2] for edge in pell_edges}) == 6
    assert len({edge[3] for edge in pell_edges}) == 6
    print(f"PASS: Pell norm N={N} has selected C={selected_C(pell, pell_rows):.12g}<2;")
    print("      its six selected pair squareclasses are all distinct.")


if __name__ == "__main__":
    main()
