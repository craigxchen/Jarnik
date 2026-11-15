# twin' (class z == c*z0 mod q, c in (Z/q)^*) = representations of N by ONE binary form of discriminant -4q^2
import math, sympy, random
from math import gcd, isqrt
from common import gpoints
random.seed(3)
split=[p for p in sympy.primerange(5,200) if p%4==1]
tot=0;bad=0
for _ in range(25):
    N=math.prod(random.sample(split[:10],random.randint(3,5))); P=gpoints(N)
    q=random.choice([l for l in sympy.primerange(3,60) if N%l])
    z0=random.choice(P)
    # direct count
    cls=set()
    for c in range(1,q):
        cls.add(((c*z0[0])%q,(c*z0[1])%q))
    direct=sum(((x%q,y%q) in cls) for x,y in P)
    # lattice L = Z z0 + q Z[i]: HNF basis from generators z0, q, q*i
    M=sympy.Matrix([[z0[0],z0[1]],[q,0],[0,q]])
    H=M.T.columnspace()  # not a basis in general; use hermite via sympy
    from sympy.matrices.normalforms import hermite_normal_form
    B=hermite_normal_form(sympy.Matrix([[z0[0],q,0],[z0[1],0,q]]))
    b1=(int(B[0,0]),int(B[1,0])); b2=(int(B[0,1]),int(B[1,1]))
    det=abs(b1[0]*b2[1]-b1[1]*b2[0]); assert det==q, (det,q)
    a=b1[0]**2+b1[1]**2; b=2*(b1[0]*b2[0]+b1[1]*b2[1]); cc=b2[0]**2+b2[1]**2
    assert b*b-4*a*cc==-4*q*q
    # count reps f(X,Y)=N : solve via points (each lattice vector z=X b1+Y b2 with |z|^2=N)
    reps=0
    for x,y in P:
        # solve x = X b1x + Y b2x, y = X b1y + Y b2y
        D=b1[0]*b2[1]-b1[1]*b2[0]
        X=(x*b2[1]-y*b2[0]); Y=(b1[0]*y-b1[1]*x)
        if X%D==0 and Y%D==0: reps+=1
    tot+=1; bad+=(reps!=direct)
print("cases",tot,"mismatches",bad)
