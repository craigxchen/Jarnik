"""Exact cube-free certificate in the constant and ramified ternary models.

Four formal marks are certified through their exact factor exponents. The
specialized squarefree polynomial whose square contains their common norm
is expanded over F3; its derivative and its gcd with the derivative are
checked in both models. Requires python-flint.
"""
from itertools import product

from flint import fmpz_mod_poly_ctx, fq_default_ctx, nmod_poly


def formal_exponents():
    fp3 = fmpz_mod_poly_ctx(3)
    fq = fq_default_ctx(modulus=fp3([1, -1, 0, 1]), var="a")
    a = fq.gen()
    elements = [u+v*a+w*a**2 for u, v, w in product(range(3), repeat=3)]
    chi = lambda x: 0 if x == 0 else 1 if x**13 == 1 else -1
    rows = [[], [], [], []]
    for x, y in product(elements, repeat=2):
        forms = (x, y, x+y, x-y)
        for row, z in zip(rows, forms):
            assert chi(z**3) == chi(z) and chi(-z) == -chi(z)
            row.append(1+chi(z))
    assert len({tuple(row) for row in rows}) == 4
    assert all(sum(row) == 729 for row in rows)
    common = sum(min(column) for column in zip(*rows))
    assert common == 92
    return common


def check(d, delta, common):
    t = nmod_poly([0, 1], 3)
    q, n, h = 27, 729, 27
    A = t**q-t
    # W consists of w0*(x+y*t^3) over F_q, where w0=A^3*eta^3, eta^2=d.
    # The two models have d=-1 and d=2t^2+1; every finite zero of d
    # lies over a zero of A.
    assert d.degree() == 0 or not (A % d)
    u = 12
    T = 1+t*A**(3*u)
    L = T.degree()
    root_height = 3*(A.degree()+delta+1)
    c = A**(3*(q-1)) * d**(3*(q-1)//2)  # w0^(q-1)
    e = c**q * A**(3*(q-1))               # P_<w0>(t^3*w0)^(q-1)
    H = T**q-c*T
    R = H**q-e*H
    assert R.degree() == n*L and R % A == 1
    assert R.derivative() == c*e*A**(3*u)
    assert R.gcd(R.derivative()) == 1
    point_degree = (n-common)*L
    max_difference_bound = ((n+h)//2-common)*L + root_height*(n-h)//2
    assert 2*max_difference_bound < point_degree
    print(f"d={d}, L={L}, deg(R)={R.degree()}, deg gcd(R,R')=0, "
          f"deg(mark)={point_degree}, deg(N)={2*point_degree}, "
          f"pair degree <= {max_difference_bound}")


if __name__ == "__main__":
    g = formal_exponents()
    print("n=729, h=27, g=", g, sep="")
    check(nmod_poly([2], 3), 0, g)
    check(nmod_poly([1, 0, 2], 3), 1, g)
