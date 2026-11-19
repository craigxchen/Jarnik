"""Exact seven-cut, monochromatic-four-cycle certificate. Standard library only."""
from fractions import Fraction
from functools import reduce
from itertools import combinations
from math import gcd, lcm

ROOTS=((0,6),(-135,-36),(-135,-51),(-135,75),(-145,81),(-117,-45),(-108,-30))
ROWS=(set(), {0,1,2,4}, {0,1,3,5}, {0,2,3,6})
CYCLE=(0,1,4,5)
RESIDUE=79658575767
PERIOD=86434978500


def mul(z,w):
    return z[0]*w[0]-z[1]*w[1], z[0]*w[1]+z[1]*w[0]


def conjugate(z):
    return z[0],-z[1]


def norm(z):
    return z[0]**2+z[1]**2


def exact_div(z,w):
    a,b=mul(z,conjugate(w))
    v=norm(w)
    assert a%v==b%v==0
    return a//v,b//v


def ggcd(z,w):
    while norm(w):
        h=norm(w)
        a,b=mul(z,conjugate(w))
        c,d=mul(((2*a+h)//(2*h),(2*b+h)//(2*h)),w)
        z,w=w,(z[0]-c,z[1]-d)
    return z


def words():
    ret=[]
    for S,T in combinations(ROWS,2):
        poly=[(1,0)]
        for j in sorted(S^T):
            re,im=ROOTS[j]
            c=(-re,-im if j in S else im)
            new=[(0,0)]*(1+len(poly))
            for k,z in enumerate(poly):
                a,b=mul(z,c)
                new[k]=(new[k][0]+a,new[k][1]+b)
                new[k+1]=(new[k+1][0]+z[0],new[k+1][1]+z[1])
            poly=new
        assert len(poly)==5 and all(b==0 for a,b in poly[1:])
        ret.append((tuple(a for a,b in poly),poly[0][1]))
    return ret


def evaluate(poly,n):
    return sum(a*n**k for k,a in enumerate(poly))


def check():
    data=words()
    masks=[sum(1<<i for i,S in enumerate(ROWS) if j in S) for j in range(7)]
    assert sorted(masks)==list(range(2,16,2))
    assert len({z for re,im in ROOTS for z in [(re,im),(re,-im)]})==14
    bs=[abs(b) for _,b in data]
    assert bs==[19966500,16110900,14871600,3855600,5094900,1239300]
    assert lcm(*(bs[i] for i in CYCLE)) == PERIOD
    divisors=[gcd(evaluate(f,RESIDUE),b) for f,b in data]
    assert [bs[i]//divisors[i] for i in CYCLE] == [2]*4
    print('quartic imaginary constants (absolute)',bs)
    print('four cycle gcds on the progression',[divisors[i] for i in CYCLE])
    n=RESIDUE
    h=[(n-re,-im) for re,im in ROOTS]
    points=[reduce(mul,(z if j in S else conjugate(z) for j,z in enumerate(h)),(1,0))
            for S in ROWS]
    den=reduce(ggcd,points)
    points=[exact_div(z,den) for z in points]
    N=norm(points[0])
    assert len(set(points))==4 and all(norm(z)==N for z in points)
    assert norm(reduce(ggcd,points)) == 1
    offsets=[]
    for z,w in combinations(points,2):
        a,b=mul(z,conjugate(w))
        half=Fraction(b,N+a)
        assert half.denominator % 2 and half.numerator % 2 == 0
        offsets.append(abs(half.numerator))
    assert offsets==[2,2,7344,560,2,2]
    tangents=[Fraction(mul(z,conjugate(points[0]))[1],N+mul(z,conjugate(points[0]))[0])
              for z in points]
    bound=2*(max(tangents)-min(tangents)) # width by the 1-Lipschitz arctangent
    print('primitive parent N',N,'digits',len(str(N)))
    print('offsets on the six lexically ordered edges',offsets)
    print('reciprocal critical arc constant exceeds',
          next(2**j for j in range(100) if (2**(j+1)*bound)**4*N>=1))


if __name__=='__main__':
    check()
