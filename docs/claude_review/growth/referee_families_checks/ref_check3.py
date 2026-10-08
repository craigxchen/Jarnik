import sys, math
sys.path.insert(0,'.')
from decimal import Decimal as D, getcontext
import ref_check as rc
getcontext().prec=120
from math import isqrt
# ---- axis six-point family: y2^2 - 2 y1^2 = -146, y1 odd, a=(y1^2-143)/2
sols=[]
for y1 in range(1,200000,2):
    r=2*y1*y1-146
    if r>0:
        t=isqrt(r)
        if t*t==r: sols.append((y1,t))
print("fundamental-range solutions (y1 odd, y1<2e5):",sols[:12])
def six(y1,y2):
    a=(y1*y1-143)//2; assert (y1*y1-143)%2==0
    pts=[(a,12),(a,-12),(a-1,y1),(a-1,-y1),(a-2,y2),(a-2,-y2)]
    n=a*a+144
    assert all(x*x+y*y==n for x,y in pts)
    return a,pts
seen=set()
for (y1,y2) in sols[:8]:
    a,pts=six(y1,y2)
    C,N,k,gn,Ps=rc.cluster_C(pts)
    print("  y1=%d y2=%d a=%d k=%d gcdnorm=%d N=%d C=%.10f"%(y1,y2,a,k,gn,N,C))
# orbit of the claimed seed under 3+2sqrt2 (parity of y1 is preserved by a single application)
y1,y2=1171,1656; out=[]
for it in range(8):
    a,pts=six(y1,y2); C,N,k,gn,Ps=rc.cluster_C(pts)
    out.append(C); print("  orbit a=%d gcdnorm=%d C-4=%.3e"%(a,gn,C-4))
    y2,y1=3*y2+4*y1,2*y2+3*y1
    assert y1%2==1 and y2*y2-2*y1*y1==-146
print("all C>4:",all(c>4 for c in out))
# ---- record clusters (independently recomputed points from the stated n)
def pts_at(n,xs):
    P=[]
    for x in xs:
        r=n-x*x; y=isqrt(r); assert y*y==r, (n,x); P+= [(x,y)] if y==0 else [(x,y),(x,-y)]
    return P
C,N,k,gn,Ps=rc.cluster_C(pts_at(32988780417125,[5743586,5743585,5743583,5743582]))
print("record8 a=5743586: k=%d gcdnorm=%d C=%.9f"%(k,gn,C))
C,N,k,gn,Ps=rc.cluster_C(pts_at(348984761607530165,[590749322-d for d in (0,1,3,4)]))
print("record8 a=590749322: k=%d gcdnorm=%d C=%.9f"%(k,gn,C))
p7=[(14650, 4555), (14629, 4622), (14555, 4850), (14554, 4853), (14491, 5038),(14453, 5146), (14443, 5174)]
assert all(x*x+y*y==235370525 for x,y in p7)
C,N,k,gn,Ps=rc.cluster_C(p7); print("record7: k=%d gcdnorm=%d C=%.9f"%(k,gn,C))
# all lattice points of circle 235370525 near that arc: list angles to check no hidden structure
n=235370525; allp=[(x,isqrt(n-x*x)) for x in range(isqrt(n)+1) if isqrt(n-x*x)**2==n-x*x]
print(" circle 235370525 has",len(allp),"points in first quadrant (x>=0,y>=0)")
p5=[(706952, 428761), (706544, 429433), (706513, 429484), (706087, 430184),(706009, 430312)]
assert all(x*x+y*y==683617125425 for x,y in p5)
C,N,k,gn,Ps=rc.cluster_C(p5); print("record5: k=%d gcdnorm=%d C=%.9f"%(k,gn,C))
