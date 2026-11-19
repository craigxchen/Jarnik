"""ref2_tri: the Vieta-reflected triple of Astra's example is a genuine triangle at a different parent."""
import sys; sys.path.insert(0,'.')
from genuine import *
from math import isqrt
from sympy import factorint
a,b,c = 1417,61,4553
for lab,off in [(a*b,1),(b*c,2),(c*a,1)]:
    s = isqrt(lab - off*off); print(lab, s*s+off*off==lab, gcd(s,off))
Np = a*b*c; print('N\'=',Np, factorint(Np))
pts = points_on(Np)
labs = {a*b, b*c, c*a}
hits = 0
for i,z0 in enumerate(pts):
    for z1 in pts:
        A = half(z1,z0,Np)
        if not A or nm(A)!=a*b or abs(A[1])!=1: continue
        for z2 in pts:
            B_ = half(z2,z0,Np)
            if not B_ or nm(B_)!=c*a or abs(B_[1])!=1: continue
            C_ = half(z2,z1,Np)
            if C_ and nm(C_)==b*c and abs(C_[1])==2:
                r = check_triple(z0,z1,z2,Np); hits += 1
                if hits==1: print('genuine reflected triangle:', r, 'arc ratio', abs(math.atan2(z2[1],z2[0])-math.atan2(z1[1],z1[0]))*math.sqrt(Np)/Np**0.25)
print('hits', hits)
