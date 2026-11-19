import sympy as sp
from sympy import I, Poly, symbols
from fractions import Fraction as F
from gcdbound import cleared_polys, t
from spec import tuple_at, ggcd, nrm
for sd in [(F(1,2),F(-2),F(7,3),F(4)),(F(-1,3),F(5,2),F(6),F(-7))]:
    zs=[Poly(z,t,domain='QQ_I') for z in cleared_polys(1,sd)]
    coeffs=[Poly(1,t,domain='QQ_I')]; h=zs[0]
    for z in zs[1:]:
        a,b,g=sp.gcdex(h,z)          # a*h+b*z=g
        coeffs=[c*a for c in coeffs]+[b]; h=g
    lc=h.LC(); h=h.quo_ground(lc); coeffs=[c.quo_ground(lc) for c in coeffs]
    assert h.degree()==0
    # denominator: common denominator of all coefficient entries (Gaussian rationals)
    den=1
    for c in coeffs:
        for co in c.all_coeffs():
            co=sp.nsimplify(co); re_,im_=sp.re(co),sp.im(co)
            den=sp.ilcm(den, sp.fraction(sp.Rational(re_))[1], sp.fraction(sp.Rational(im_))[1])
    chk=sp.expand(sum((c*z).as_expr() for c,z in zip(coeffs,zs)))
    print([str(x) for x in sd],'sum a_i z_i =',chk,' rational denominator c =',sp.factorint(den))
    bad=[T for T in range(-300,301) if (lambda z:(lambda g: g!=(0,0) and (den*den)%nrm(g)!=0)(ggcd(ggcd(ggcd(ggcd(z[0],z[1]),z[2]),z[3]),z[4])))(tuple_at(1,sd,T))]
    print('  T with N(gcd) not dividing c^2:',bad)
