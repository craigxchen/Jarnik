# Twin: q/Q_u | F_u on affine t-flats for bits = Boolean functions of <= d affine forms (incl. p=2 powers).
# Circle: exact Gaussian identity prod_x z_x^{chi_u(x)} = (G_u/conj G_u)^{m0}, and Im G_u != 0 whenever G_u/conj G_u != 1.
import random
from itertools import product
from sympy.ntheory.modular import crt
from sympy import primerange
random.seed(1)
def dot(u,x): return bin(u&x).count('1')&1
r=6; n=1<<r
def randbit(d):
    ls=[random.randrange(1,n) for _ in range(d)]; cs=[random.randrange(2) for _ in range(d)]
    g=random.randrange(1<<(1<<d))
    return lambda x:(g>>sum(((dot(l,x)^c)<<k) for k,(l,c) in enumerate(zip(ls,cs))))&1
def rank2(vs):
    rows=[]
    for v in vs:
        for rr in rows: v=min(v,v^rr)
        if v: rows.append(v)
    return len(rows)
ppows=[8,9,25,7,11,13,17,19,23,29,31,37,41,43,4*0+47,16*0+53]
d=2; t=4
bits=[randbit(d) for _ in ppows]
cls=[(random.randrange(pp),random.randrange(pp)) for pp in ppows]
q=1
for pp in ppows: q*=pp
y={x:int(crt(ppows,[cls[i][bits[i](x)] for i in range(len(ppows))])[0]) + q*random.randrange(5) for x in range(n)}
bad=0;tests=0
for trial in range(300):
    while True:
        x0=random.randrange(n); M=[random.randrange(1,n) for _ in range(t)]
        if rank2(M)==t: break
    pts=[x0^ (M[0]*0) for _ in range(1)]
    flat=[]
    for z in range(1<<t):
        x=x0
        for k in range(t):
            if (z>>k)&1: x^=M[k]
        flat.append(x)
    for u in range(1,1<<t):
        F=sum((-1)**dot(u,z)*y[flat[z]] for z in range(1<<t))
        Qu=1
        for i,pp in enumerate(ppows):
            S=sum((-1)**dot(u,z)*bits[i](flat[z]) for z in range(1<<t))
            if S: Qu*=pp
        tests+=1
        if F % (q//Qu): bad+=1
    # support dimension <= d per bit
    for i in range(len(ppows)):
        sup=[u for u in range(1,1<<t) if sum((-1)**dot(u,z)*bits[i](flat[z]) for z in range(1<<t))]
        assert rank2(sup)<=d
print("twin divisibility tests",tests,"failures",bad)
# circle
ps=[p for p in primerange(5,200) if p%4==1][:10]
def gp(p):
    for a in range(1,p):
        b2=p-a*a
        if b2>0:
            b=int(round(b2**0.5))
            if b*b==b2: return (a,b)
def mul(z,w): return (z[0]*w[0]-z[1]*w[1], z[0]*w[1]+z[1]*w[0])
def conj(z): return (z[0],-z[1])
def pw(z,k):
    o=(1,0)
    for _ in range(k): o=mul(o,z)
    return o
pis=[gp(p) for p in ps]; E=[random.randrange(1,4) for _ in ps]
thr=[[randbit(d) for j in range(E[i])] for i in range(len(ps))]
def a(x,i):  # exponent = number of thresholds; make thresholds nested by sorting? use count of bits as exponent (unary need nested)
    return sum(thr[i][j](x) for j in range(E[i]))
# NOTE: unary thresholds must be nested; general exponent = sum of 1_{a>=j}. Using a = sum of arbitrary bits gives a valid
# exponent but bits 1_{a>=j} are then not the given functions; for the identity check only a_xp = sum_j b_j is used.
eta=(0,1)
def z(x):
    o=eta
    for i,pi in enumerate(pis):
        o=mul(o,mul(pw(pi,a(x,i)),pw(conj(pi),E[i]-a(x,i))))
    return o
m0=1<<(t-d); bad=0; tests=0; imz=0
for trial in range(100):
    while True:
        x0=random.randrange(n); M=[random.randrange(1,n) for _ in range(t)]
        if rank2(M)==t: break
    flat=[]
    for zz in range(1<<t):
        x=x0
        for k in range(t):
            if (zz>>k)&1: x^=M[k]
        flat.append(x)
    for u in range(1,1<<t):
        num=(1,0);den=(1,0)
        for zz in range(1<<t):
            if dot(u,zz): den=mul(den,z(flat[zz]))
            else: num=mul(num,z(flat[zz]))
        G=(1,0)
        for i,pi in enumerate(pis):
            s=sum((-1)**dot(u,zz)*thr[i][j](flat[zz]) for zz in range(1<<t) for j in range(E[i]))
            assert s%m0==0
            b=s//m0
            G=mul(G,pw(pi,b) if b>0 else pw(conj(pi),-b))
        # check num/den == (G/conjG)^{m0}: num*conj(G)^m0 == den*G^m0
        tests+=1
        if mul(num,pw(conj(G),m0))!=mul(den,pw(G,m0)): bad+=1
        if G[1]==0 and G!=(1,0) : imz+=1
print("circle identity tests",tests,"failures",bad,"real nontrivial G",imz)
