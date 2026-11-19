"""Build my own integer specialisation of the V.3 four-block star (own parameter, small),
for use as base model in the V.2 full-isotropy construction."""
from fractions import Fraction as F
from itertools import combinations
from math import gcd, prod
from sympy import factorint
roots=(1,-9,-7,2,4); u=[r*r for r in roots]
w=roots[1:]
s1=sum(w); s2=sum(a*b for a,b in combinations(w,2)); s3=sum(a*b*c for a,b,c in combinations(w,3)); s4=prod(w)
assert s3==s1*s2 and s1*s4!=0
edges=list(combinations(range(5),2))
def model(D,P):
    t=F(P,D)
    Gr=t**4-s2*t**2+s4; Gi=-s1*t**3+s3*t
    Ar=-t/s1; Ai=F(1)
    cols=[]
    for ui in u:
        m=F(ui-1,ui)*(ui+t*t)
        re=D**4*(Gr+m*Ar); im=D**4*(Gi+m*Ai)
        if re.denominator!=1 or im.denominator!=1: return None
        cols.append((int(re),int(im)))
    blocks=[P*P+(D*roots[i]*roots[j])**2 for i,j in edges]
    cof=[]
    for (i,j),J in zip(edges,blocks):
        a,b=cols[i],cols[j]; d=a[0]*b[1]-a[1]*b[0]
        if d%J: return None
        cof.append(d//J)
    return cols,blocks,cof
D=6350400
for P in range(700*D+1, 700*D+4000000, 2):
    if gcd(P,D)!=1: continue
    r=model(D,P)
    if r is None: continue
    cols,blocks,cof=r
    if any(gcd(*c)!=1 for c in cols): continue
    Np=prod(blocks)
    if any(gcd(a,b)!=1 for a,b in combinations(blocks,2)): continue
    if any(gcd(abs(c),Np)!=1 for c in cof): continue
    fs=[factorint(J) for J in blocks]
    if any(e>1 or p%4!=1 for f in fs for p,e in f.items()): continue
    print("P=",P,"D=",D); print("cof",cof); print("blocks",blocks); print("factors",fs)
    break
