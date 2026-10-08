"""Independent exact-domain and numerical-arc audit for pair rotations."""

from itertools import combinations
from math import atan2, pi, sqrt
from random import Random

from check_quartet_matching_gcd_cut_budget import ggcd, mul, norm, sub


def quotient_if_integral(a, b):
    d = norm(b)
    num = mul(a, (b[0], -b[1]))
    if num[0] % d or num[1] % d:
        return None
    return num[0] // d, num[1] // d


def audit(points, angles, circles):
    N = norm(points[0])
    R = sqrt(N)
    L = R * (angles[-1] - angles[0])
    C = L / sqrt(R)
    for ia in range(len(points)):
        for ib in range(len(points)):
            if ia == ib:
                continue
            a, b = points[ia], points[ib]
            common = ggcd(a, b)
            A, B = quotient_if_integral(b, common), quotient_if_integral(a, common)
            H = norm(B)
            assert norm(A) == H and norm(ggcd(A, B)) == 1 and N % H == 0
            chord2, reduced_chord2 = norm(sub(b, a)), norm(sub(A, B))
            assert reduced_chord2 >= 2 and reduced_chord2 % 2 == 0
            assert H * chord2 == N * reduced_chord2
            divided = []
            for z in points:
                u = quotient_if_integral(z, B)
                rotated = quotient_if_integral(mul(A, z), B)
                assert (u is None) == (rotated is None)
                if u is not None:
                    assert norm(u) == N // H and norm(rotated) == N
                    divided.append(u)
            assert a in [z for z in points if quotient_if_integral(z, B) is not None]
            assert all(norm(sub(u, v)) >= 2 and norm(sub(u, v)) % 2 == 0
                       for u, v in combinations(divided, 2))
            full_domain = [z for z in circles[N] if quotient_if_integral(z, B) is not None]
            assert len(full_domain) == len(circles[N // H])
            count = len(divided)
            gap = R * abs(angles[ib] - angles[ia])
            assert count <= 1 + L / sqrt(2 * H) + 1e-10
            assert count <= 1 + L * sqrt(chord2) / (2 * R) + 1e-10
            assert count < 1 + C * C * gap / (2 * L) + 1e-10
            if C <= sqrt(2):
                assert count == 1
    return len(points) * (len(points) - 1)


def audit_codes():
    for m in (4, 6, 8, 10):
        cuts = [set(cut) for cut in combinations(range(m), m // 2)]
        for a in range(m):
            for b in range(m):
                if a == b:
                    continue
                rotation = [z for z in range(m) if all(
                    (z in cut) == (a in cut) for cut in cuts if (a in cut) != (b in cut))]
                reflection = [z for z in range(m) if all(
                    (z in cut) == (a in cut) for cut in cuts if (a in cut) == (b in cut))]
                assert rotation == [a]
                assert set(reflection) == {a, b}


def main():
    circles = {}
    for x in range(-31, 32):
        for y in range(-31, 32):
            N = x * x + y * y
            if 0 < N <= 1000:
                circles.setdefault(N, []).append((x, y))
    ordered = {}
    for N, points in circles.items():
        rows = sorted((atan2(y, x), (x, y)) for x, y in points)
        ordered[N] = rows + [(angle + 2 * pi, z) for angle, z in rows]
    rng, count = Random(20260917), 0
    keys = list(circles)
    for _ in range(1000):
        N = rng.choice(keys)
        start = rng.randrange(len(circles[N]))
        length = rng.randrange(2, min(8, len(circles[N])) + 1)
        window = ordered[N][start:start + length]
        count += audit([z for angle, z in window], [angle for angle, z in window], circles)
    audit_codes()
    print(f"PASS: {count:,} ordered anchor-pair domains in 1,000 circle arcs;")
    print("full-circle domain counts, exact arithmetic, capacities, balanced codes through m=10.")


if __name__ == "__main__":
    main()
