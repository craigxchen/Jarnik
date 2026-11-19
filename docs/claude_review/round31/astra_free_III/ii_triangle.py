"""Exact common-parent test for the Heron equation and its Vieta move."""
from math import gcd, isqrt, prod


def main():
    axes = (1417, 20381, 4553)
    a, b, c = axes
    u, v, w = 2, 1, 1
    parent = prod(axes)
    for norm, offset in ((b*c, u), (c*a, v), (a*b, w)):
        real = isqrt(norm)
        assert real**2 + offset**2 == norm and gcd(real, offset) == 1
        assert norm**2 >= 256*offset**4*parent and parent % norm == 0
    A, B, C = a*u*u, b*v*v, c*w*w
    U = (u*v*w)**2
    assert (A+B-C)**2 == 4*A*B-4*U
    s = isqrt(a*b)
    assert c*w*w in (A+B-2*s*u*v, A+B+2*s*u*v)
    sign = 1 if c*w*w == A+B+2*s*u*v else -1
    assert w*w*(a*c-v*v) == (a*u+sign*s*v)**2
    assert w*w*(b*c-u*u) == (s*u+sign*b*v)**2
    assert gcd(gcd(a,b),c) == 1
    reflected = 2*(A+C)-B
    assert reflected == 61 and parent % reflected == 40
    assert (A+reflected-C)**2 == 4*A*reflected-4*U
    print(f"parent={parent}, squared lengths={(A,B,C)}, U={U}")
    print(f"reflected largest length={reflected}; parent remainder={parent % reflected}")


if __name__ == '__main__':
    main()
