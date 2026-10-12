"""Proposition 6.3: Cayley coordinates of a normalised configuration and reduced chords.

For random circles N (products of split primes) and all pairs of lattice points z, w:
  t = i(rho+1)/(rho-1) with rho = z/w is rational; denominator(t) vs beta = (z - w)/gcd(z, w);
  the table: denominator in {|beta|/2, |beta|/sqrt2}; every odd rational prime dividing beta divides it.
Also checks the addition formula T_xy = (t_x t_y + 1)/(t_y - t_x) for triples.  Exact (Fractions).
"""
import random
from fractions import Fraction
import sympy
from g8common import gmul, gconj, gnorm, gdivexact, ggcd, gauss_prime


def cayley(z, w):
    """t = i (rho + 1)/(rho - 1), rho = z/w; computed as i (z + w)/(z - w)."""
    num = gmul((0, 1), (z[0] + w[0], z[1] + w[1]))
    den = (z[0] - w[0], z[1] - w[1])
    q = gmul(num, gconj(den))
    n = gnorm(den)
    assert q[1] == 0, "t not real"
    return Fraction(q[0], n)


def points(primes):
    pis = [gauss_prime(p) for p in primes]
    pts = [(1, 0)]
    for pi in pis:
        pts = [gmul(z, f) for z in pts for f in (pi, gconj(pi))]
    out = []
    for z in pts:
        u = z
        for _ in range(4):
            out.append(u)
            u = (-u[1], u[0])
    return out


def main():
    rng = random.Random(7)
    split = [p for p in sympy.primerange(5, 2000) if p % 4 == 1]
    npairs = 0
    ntriples = 0
    for trial in range(40):
        primes = rng.sample(split, rng.randint(2, 5))
        pts = points(primes)
        for _ in range(60):
            z, w = rng.sample(pts, 2)
            g = ggcd(z, w)
            beta = gdivexact((z[0] - w[0], z[1] - w[1]), g)
            t = cayley(z, w)
            den = t.denominator
            b2 = gnorm(beta)
            assert 4 * den * den == b2 or 2 * den * den == b2, (z, w, t, beta)
            for l in sympy.primefactors(b2):
                if l % 2 == 1 and gdivexact(beta, (l, 0)) is not None:
                    assert den % l == 0
            npairs += 1
        for _ in range(20):
            z1, zx, zy = rng.sample(pts, 3)
            tx, ty = cayley(zx, z1), cayley(zy, z1)
            T = cayley(zx, zy)
            assert T == (tx * ty + 1) / (ty - tx)
            ntriples += 1
    print(f"pairs checked: {npairs} (denominator = |beta|/2 or |beta|/sqrt2; odd primes of beta divide it): PASS")
    print(f"triples checked: {ntriples} (T_xy = (t_x t_y + 1)/(t_y - t_x)): PASS")


if __name__ == '__main__':
    main()
