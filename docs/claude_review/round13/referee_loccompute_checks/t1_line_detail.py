import sys, itertools, math
from rlib import *
L=int(sys.argv[1]); D=int(sys.argv[2]); K=int(sys.argv[3])
P3=[p for p in itertools.product(range(-L,L+1),repeat=3) if sum(p)%2==0]
ys=sp.symbols('a b c')
JP=J_exact(P3,D)
def prim(v):
    g=0
    for x in v: g=math.gcd(g,abs(x))
    return g==1
dirs=[v for v in itertools.product(range(-K,K+1),repeat=3) if prim(v) and next(x for x in v if x)>0]
W=[]
for v in dirs:
    vals={sum(a*b for a,b in zip(v,p)) for p in P3}
    if len(vals)<=D: W.append(sp.Poly(sp.expand(sum(a*y for a,y in zip(v,ys))**len(vals)),*ys))
base=ideal_dim(W,ys,D)[0]
print("width",base,"J",len(JP[D]))
G=list(W)
for n in dirs:
    ns=sp.Matrix([list(n)]).nullspace()
    B=[]
    for v in ns:
        den=1
        for x in v: den=sp.ilcm(den,sp.fraction(x)[1])
        B.append([int(x*den) for x in v])
    img=sorted({tuple(sum(a*b for a,b in zip(u,p)) for u in B) for p in P3})
    Jimg=J_exact(img,D)
    new=[]
    for DD in range(D+1):
        if Jimg[DD]: new+=pull_forms(Jimg[DD], monos(2,DD), B, ys)
    r=ideal_dim(W+new,ys,D)[0]
    if r>base:
        print("kernel n",n,"basis",B,"|img|",len(img),"Jimg dims",[len(Jimg[i]) for i in range(D+1)],"alone gives",r-base)
    G+=new
print("all",ideal_dim(G,ys,D)[0])
