# Independent check of B^4=70: count lattice points of n*R_5, R_5={x in Z^5: sum x=0, sum |x| <= 2}
from fractions import Fraction
from math import comb
import itertools
def count(n):
    # number of x in Z^5, sum x = 0, sum |x| <= 2n
    # enumerate x1..x4, x5 = -sum
    c=0
    r=2*n
    for x1 in range(-r,r+1):
        for x2 in range(-r+abs(x1), r-abs(x1)+1):
            s2=abs(x1)+abs(x2)
            for x3 in range(-(r-s2), r-s2+1):
                s3=s2+abs(x3)
                for x4 in range(-(r-s3), r-s3+1):
                    s4=s3+abs(x4)
                    x5=-(x1+x2+x3+x4)
                    if s4+abs(x5)<=r: c+=1
    return c
vals=[count(n) for n in range(0,7)]
print(vals)
# fit degree-4 polynomial via finite differences
import sympy as sp
n=sp.symbols('n')
poly=sp.interpolate(list(zip(range(7),vals)),n)
print(sp.expand(poly), 'B^4 =', sp.Poly(poly,n).LC()*24)
# Astra's formula
print(sum(comb(5,k)*comb(3,k-1) for k in range(1,5)))
# Astra's off-hyperplane count for n
def astra(nn):
    return sum(comb(5,k)*sum(comb(s-1,k-1)*comb(s-1,4-k) for s in range(1,nn+1)) for k in range(1,5))
print([astra(m) for m in range(7)])
