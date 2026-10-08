"""Literal circle and angular fixtures for the virtual-normal constraints."""

from itertools import combinations
from math import atan2, cos, gcd, hypot, isqrt, pi, sin, sqrt, tan


def circle_points(n):
    points = set()
    for x in range(-isqrt(n), isqrt(n)+1):
        y = isqrt(n-x*x)
        if x*x+y*y == n:
            points.update(((x, y), (x, -y)))
    return sorted(points, key=lambda z: atan2(z[1], z[0]))


def main():
    pairs = windows = short_anchors = projections = 0
    for n in (5, 13, 25, 50, 65, 100, 125):
        radius = sqrt(n)
        points = circle_points(n)
        for z, w in combinations(points, 2):
            cross, dot = z[0]*w[1]-z[1]*w[0], z[0]*w[0]+z[1]*w[1]
            if cross == 0:
                continue
            gz, gw = gcd(*z), gcd(*w)
            assert cross % (gz*gw) == 0 and abs(cross) >= gz*gw
            angle = atan2(abs(cross), dot)
            assert (n/(gz*gw))*angle >= 1-1e-12
            pairs += 1
        angles = [(atan2(z[1], z[0]), z) for z in points]
        cyclic = angles + [(angle+2*pi, z) for angle, z in angles]
        for angle, z in angles:
            cohort = [w for beta, w in cyclic
                      if angle <= beta <= angle+2/sqrt(radius)]
            # C=2: rho<R^(1/4)/sqrt(2) iff 16*(N/g^2)^4<N.
            short = [w for w in cohort if 16*(n//gcd(*w)**2)**4 < n]
            assert len(short) <= 1
            short_anchors += len(short)
            windows += 1
    assert short_anchors > 0

    for proposed_constant in (1, 2, 5, 10):
        radius = (8*proposed_constant)**4
        epsilon = 1/sqrt(radius)
        assert epsilon <= 0.1 and tan(3*epsilon) <= 4*epsilon
        assert sqrt(radius)/4 > proposed_constant*radius**0.25
        length = 2*sqrt(radius)
        bound = length*abs(sin(2*epsilon))+length*length/(4*radius)
        assert bound <= 5

    for radius in (2, 5, 100, 10000):
        for length in (0.1*sqrt(radius), sqrt(radius), 2*sqrt(radius)):
            for delta in (0, 0.1, 0.7, pi/2, 2.4):
                for v in ((1, 0), (2, 1), (3, 4)):
                    rho, a = hypot(*v), length/(2*radius)
                    values = [rho*radius*cos(delta-a+2*a*j/100)
                              for j in range(101)]
                    bound = rho*(length*abs(sin(delta))+length*length/(4*radius))
                    assert max(values)-min(values) <= bound+1e-9
                    projections += 1
    print(f'PASS: {pairs} actual determinant pairs; {windows} endpoint windows '
          f'({short_anchors} short-anchor incidences); 4 angular obstructions; '
          f'{projections} projection-width fixtures.')


if __name__ == '__main__':
    main()
