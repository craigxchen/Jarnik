"""Search rational k=4 polynomial (degree-1 block) families: roots rho=i, a=r12, b=r13, c=r23 in Q(i) with
Phi(i,x,y)=0 for all pairs among {a,b,c}, where
Phi(i,(u2,v2),(u3,v3)) = v2 v3 (u2-u3)^2 + v2 u2^2 + v3 u3^2 + (1+v2)(1+v3)(v2+v3)(1+v2+v3)."""
from fractions import Fraction as F
from math import isqrt
import itertools, sys
def phi(a,b):
    (u2,v2),(u3,v3)=a,b
    return v2*v3*(u2-u3)**2 + v2*u2**2 + v3*u3**2 + (1+v2)*(1+v3)*(v2+v3)*(1+v2+v3)
def ratsqrt(x):
    if x < 0: return None
    n, d = x.numerator, x.denominator
    a, b = isqrt(n), isqrt(d)
    return F(a,b) if a*a==n and b*b==d else None
def neighbours(a, vgrid):
    u2,v2 = a; out=[]
    for v3 in vgrid:
        A = v2*v3 + v3; B = -2*v2*v3*u2; C = v2*v3*u2**2 + v2*u2**2 + (1+v2)*(1+v3)*(v2+v3)*(1+v2+v3)
        if A == 0: continue
        D = B*B-4*A*C; s = ratsqrt(D)
        if s is None: continue
        for u3 in {(-B+s)/(2*A), (-B-s)/(2*A)}:
            b=(u3,v3); assert phi(a,b)==0; out.append(b)
    return out
if __name__ == "__main__":
    H = int(sys.argv[1]) if len(sys.argv)>1 else 6
    vals = sorted({F(p,q) for q in range(1,H+1) for p in range(-3*H,3*H+1)} - {F(0), F(-1)})
    ugrid = sorted({F(p,q) for q in range(1,H+1) for p in range(-3*H,3*H+1)})
    found = []
    for a in itertools.product(ugrid, vals):
        if a[1] == 0: continue
        N = neighbours(a, vals)
        N = [b for b in N if b != a and b != (F(0),F(1))]
        for b, c in itertools.combinations(N, 2):
            if b != c and phi(b,c)==0:
                found.append((a,b,c))
        if len(found) >= 20: break
    for t in found[:20]: print([ (str(x[0]),str(x[1])) for x in t ])
    print("found", len(found))
