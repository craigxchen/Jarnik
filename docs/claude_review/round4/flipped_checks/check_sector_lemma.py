"""Exhaustive check of the sector lemma (Lemma 6.1 of flipped.md) on random sectors.

A closed sector {z : |arg z - t| <= eps, 0 < |z|^2 <= X} with 0 < 2 eps < pi contains at most
1 + 2 eps X primitive lattice points, and any two of them u != v satisfy N(u) N(v) >= 1/sin(2 eps)^2.
Angles are compared exactly: |arg z - t| <= eps is tested with exact integer cross products against
the two boundary rays given by rational approximations (we use rational boundary directions, so the
test is exact for the sector they define, whose half-width eps is computed from them).
"""
import math, random
from fractions import Fraction


def main():
    rng = random.Random(7)
    tested = 0
    for trial in range(400):
        X = rng.choice([200, 1000, 5000, 20000])
        # boundary rays r1, r2: integer vectors; sector = cone between them (angle < pi/2)
        while True:
            r1 = (rng.randint(-300, 300), rng.randint(-300, 300))
            ang = rng.uniform(1e-4, 0.3)
            th = math.atan2(r1[1], r1[0]) + ang
            r2 = (round(1000 * math.cos(th)), round(1000 * math.sin(th)))
            if r1 != (0, 0) and r1[0] * r2[1] - r1[1] * r2[0] > 0:
                break
        a1 = math.atan2(r1[1], r1[0]); a2 = math.atan2(r2[1], r2[0])
        two_eps = (a2 - a1) % (2 * math.pi)
        assert 0 < two_eps < math.pi / 2
        pts = []
        R = math.isqrt(X)
        for a in range(-R, R + 1):
            for b in range(-R, R + 1):
                if (a, b) == (0, 0) or a * a + b * b > X or math.gcd(a, b) != 1:
                    continue
                # inside cone: cross(r1, z) >= 0 and cross(z, r2) >= 0 (exact)
                if r1[0] * b - r1[1] * a >= 0 and a * r2[1] - b * r2[0] >= 0:
                    pts.append((a, b))
        bound = 1 + two_eps * X          # 1 + 2 eps X with 2 eps = two_eps
        assert len(pts) <= bound + 1e-9, (len(pts), bound)
        s2 = math.sin(two_eps) ** 2
        for i in range(len(pts)):
            for j in range(i + 1, len(pts)):
                u, v = pts[i], pts[j]
                nu, nv = u[0] ** 2 + u[1] ** 2, v[0] ** 2 + v[1] ** 2
                det = u[0] * v[1] - u[1] * v[0]
                assert det != 0
                assert nu * nv * s2 >= 1 - 1e-12
        tested += 1
    print("sectors tested:", tested, "PASS")


if __name__ == "__main__":
    main()
