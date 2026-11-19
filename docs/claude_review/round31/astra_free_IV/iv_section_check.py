"""Exact checks for the surface inequality in Part IV (no symbolic packages)."""

from fractions import Fraction as Q
from math import comb


a = Q(2, 3)


def contribution(k, previous, new):
    """Change of 3e_orb - (K+boundary)^2, including old pair contacts."""
    return (3 * (1 - 2*new + new*previous)
            + 3*a*k*(new - previous)
            + (1 + new - previous - a*k)**2
            - 2*a*a*comb(k, 2))


def savings(k):
    return max(0, comb(k, 2) + 2*k - 9)


def main():
    for k in range(3, 257):
        assert contribution(k, 0, 0) == 4 - 8*Q(k, 9)
    assert contribution(2, 0, a) == Q(8, 9)
    assert contribution(2, a, a) == -Q(4, 9)
    for k in range(1, 128):
        for j in range(1, 128):
            assert savings(k+j)-savings(k)-savings(j) <= 2*k*j
    assert savings(6)-2*savings(3) == 2*3*3

    for m in range(4, 161):
        for s in range(m+1):
            upper_self = s-Q(m, 4)
            raw = (Q(8,9)*upper_self + Q(4,9)*(Q(m-1,2)*upper_self-comb(s,2))
                   + 4 - Q(8,9)*s)
            assert raw == 4-Q(m,6)-Q(2,9)*(s-Q(m,2))**2
    assert Q(3*(48-24),4*48*(48+3)) == Q(1,136)
    assert Q(36-18,2*36*35) == Q(1,140)

    # A real arrangement from Part III. Start with the graphs of
    # f_i = u_i - t^2/u_i, with the products u_i*u_j all distinct. Every
    # finite pair intersects transversely twice. At infinity all m sections
    # collide twice. Blow up these last two points with new boundary weight 0.
    for m in (5,6,24):
        j = 2*comb(m,2)
        total_self = 4*m
        e = 4+2-2*a*m+a*a*j
        square = (8-4*a*m+(a*a-2*a)*total_self
                  + 2*a*a*4*comb(m,2)-2*(1-a*m)**2)
        # Four generic half-fibers cancel the constant 4-2*a*m in the ledger.
        slack = 3*e-square + 4*(-1+a*m/2)
        assert slack == Q(8,9)*total_self+a*a*j+2*contribution(m,0,0)
    print('Exact local blowup, contact-loss, and global coefficients verified.')
    print('Two complete infinity clusters and all finite intersections reproduce the direct Chern counts.')


if __name__ == '__main__':
    main()
