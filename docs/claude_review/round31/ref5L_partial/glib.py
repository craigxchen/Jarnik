from math import gcd, prod
from fractions import Fraction
import sympy as sp
def gdivmod(a,b):
    # a,b Gaussian ints as tuples; nearest-rounding quotient
    ar,ai=a; br,bi=b; n=br*br+bi*bi
    xr=ar*br+ai*bi; xi=ai*br-ar*bi
    def rnd(z): return (2*z+n)//(2*n)
    qr,qi=rnd(xr),rnd(xi)
    rr=ar-(qr*br-qi*bi); ri=ai-(qr*bi+qi*br)
    return (qr,qi),(rr,ri)
def ggcd(a,b):
    while b!=(0,0):
        _,r=gdivmod(a,b); a,b=b,r
    return a
def gmul(a,b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def gconj(a): return (a[0],-a[1])
def gexactdiv(a,b):
    q,r=gdivmod(a,b); assert r==(0,0); return q
def nrm(a): return a[0]**2+a[1]**2
def det(a,b): return a[0]*b[1]-a[1]*b[0]
def dot(a,b): return a[0]*b[0]+a[1]*b[1]
def chord2(a,b):
    # squared chord between doubled directions a^2/|a|^2 and b^2/|b|^2 : |a^2 bbar^2 - ...| ; 4 sin^2 cos^2... exact: |e^{2ia}-e^{2ib}|^2 = 4 sin^2(a-b)
    # sin^2(a-b) = det^2/(|a|^2|b|^2)
    return Fraction(4*det(a,b)**2, nrm(a)*nrm(b))
