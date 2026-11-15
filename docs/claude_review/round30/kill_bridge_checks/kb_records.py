# Independent: code of the 5.754 axis cluster (sympy Gaussian factorisation); factor R^2 of the 74-digit example.
import math, itertools, sympy as sp
from fractions import Fraction as Fr
from kb_gauss import *
def axis(a,y0,ds):
    P=[]
    for d in ds:
        y2=y0*y0+2*d*a-d*d; yy=math.isqrt(y2); assert yy*yy==y2
        P+=[(a-d,yy),(a-d,-yy)]
    return P
def vpi(z,pi):
    n=pi[0]**2+pi[1]**2; v=0
    while True:
        x=z[0]*pi[0]+z[1]*pi[1]; yy=z[1]*pi[0]-z[0]*pi[1]
        if x%n or yy%n: return v
        z=(x//n,yy//n); v+=1
def code(pts):
    n=pts[0][0]**2+pts[0][1]**2
    assert all(x*x+y*y==n for x,y in pts) and len(set(pts))==len(pts)
    f=sp.factorint(n); cols=[]
    for p,e in f.items():
        assert p%4==1
        x=next(x for x in range(1,math.isqrt(p)+1) if math.isqrt(p-x*x)**2==p-x*x); pi=(x,math.isqrt(p-x*x))
        ex=[vpi(z,pi) for z in pts]
        for tau in range(1,e+1): cols.append((p,tuple(int(v>=tau) for v in ex)))
    return n,f,cols
pts=axis(590749322,12809,(0,1,3,4))
n,f,cols=code(pts)
R=math.sqrt(n); ang=sorted(math.atan2(y,x) for x,y in pts)
print("n",n,f,"C=",R*(ang[-1]-ang[0])/math.sqrt(R))
cuts={}
for p,c in cols:
    k=min(c,tuple(1-x for x in c)); cuts[k]=cuts.get(k,[])+[p]
pairs=[(0,1),(2,3),(4,5),(6,7)]  # conjugate pairs (a-d, +-y)
print("#cuts",len(cuts))
trans=set()
for s in itertools.product([0,1],repeat=4):
    v=[0]*8
    for (i,j),b in zip(pairs,s): v[i]=b; v[j]=1-b
    v=tuple(v); trans.add(min(v,tuple(1-x for x in v)))
print("cuts == all 8 transversal cuts:", set(cuts)==trans)
M=sp.Matrix([[c[j]-c[0] for _,c in cols] for j in range(8)])
print("rank",sp.Matrix([[cols[k][1][j]-cols[k][1][0] for k in range(len(cols))] for j in range(1,8)]).rank())
# 74-digit example
y=[Fr(-9,4),Fr(-2,3),Fr(1,4),Fr(8,3)]
C,ws,R2,_=C_of(roots_cube(y),Fr(13999,2))
fac=sp.factorint(R2)
print("R2 digits",len(str(R2)),"#primes",len(fac),"squarefree",all(e==1 for e in fac.values()),"all 1 mod 4",all(p%4==1 for p in fac))
print("C",C)
# compare with lane's example.txt first point
print(ws[:2])
