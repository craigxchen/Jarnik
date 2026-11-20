# Independent exact check of the IX.5 Schottky construction (referee code; no import of Astra's scripts)
import random, math, sys
from sympy import primerange, sqrt_mod
from sympy.ntheory.modular import crt
from fractions import Fraction
random.seed(int(sys.argv[1]) if len(sys.argv)>1 else 1)
pairs=[(a,b) for a in range(-5,6) for b in range(-5,6) if (a,b)!=(0,0)]
ps=[p for p in primerange(51,20000) if p%4==1][:len(pairs)]
ell={pr:p for pr,p in zip(pairs,ps)}
res=[];mods=[]
for (a,b),l in ell.items():
    u,v=a,7*a+b; assert 0<abs(v)<=40
    s=sqrt_mod(l-1,l)        # sqrt(-1) mod l
    # u+Lv = s v  => L = (s v - u)/v
    L=((s*v-u)*pow(v,-1,l))%l
    res.append(L);mods.append(l)
L=int(crt(mods,res)[0])
def B(a,b): return (a+L*(7*a+b), 7*a+b)
for (a,b),l in ell.items():
    x,y=B(a,b); assert (x*x+y*y)%l==0
Q=1
for l in mods: Q*=l
n1=math.isqrt(sum(t*t for t in B(1,0)))+1; n2=math.isqrt(sum(t*t for t in B(0,1)))
ratio=n1//n2+2
D=Q*((102+ratio)//Q+1)
print("digits L",len(str(L)),"digits Q",len(str(Q)),"digits D",len(str(D)))
def det(u,v): return u[0]*v[1]-u[1]*v[0]
def nrm(u): return u[0]**2+u[1]**2
def word(N):
    w=[random.randint(1,3)]
    while len(w)<N:
        c=random.randint(1,3)
        if c!=w[-1]: w.append(c)
    return w
ok=True
for trial in range(20):
    w=word(14)
    P=[B(0,1),B(1,w[0]*D)]
    for n in range(1,len(w)):
        c=(w[n]-w[n-1])*D
        P.append((c*P[-1][0]+P[-2][0],c*P[-1][1]+P[-2][1]))
    for n in range(len(P)-1):
        assert abs(det(P[n],P[n+1]))==1
        assert nrm(P[n+1])>=10**4*nrm(P[n])
        if n>=1: assert nrm(P[n+1])<=(2*D+1)**2*nrm(P[n])
    for n in range(len(P)):
        ref=B(0,1) if n%2==0 else B(1,0)
        assert (P[n][0]-ref[0])%Q==0 and (P[n][1]-ref[1])%Q==0
    # small combinations A P_k + B' P_{k-1}: norm divisible by a chosen ell
    for k in range(1,len(P)):
        for A in range(-5,6):
            for Bp in range(-5,6):
                if (A,Bp)==(0,0): continue
                v=(A*P[k][0]+Bp*P[k-1][0],A*P[k][1]+Bp*P[k-1][1])
                N=nrm(v)
                if not any(N%l==0 for l in mods): ok=False; print("FAIL combo",k,A,Bp)
    # independent: CF of slope of last vector, all w=(a,b) with t^2|sin(arg w - phi)|<1, t>=h_1
    U=P[-1]; hU2=nrm(U)
    # candidates: convergents and intermediate fractions of slope U[1]/U[0]
    x,y=U; cands=set()
    # CF of y/x
    aa=[];X,Y=x,y
    if X<0: X,Y=-X,-Y
    num,den=Y,X
    while den:
        q=num//den; aa.append(q); num,den=den,num-q*den
    p0,q0,p1,q1=1,0,aa[0],1
    conv=[(q0,p0),(q1,p1)]
    for q in aa[1:]:
        for j in set(range(0,min(q,8)+1))|set(range(max(0,q-8),q+1)):
            cands.add((j*q1+q0,j*p1+p0))
        p0,q0,p1,q1=p1,q1,q*p1+p0,q*q1+q0; conv.append((q1,p1))
    for (qq,pp) in conv: cands.add((qq,pp))
    h1sq=nrm(P[1]); cnt=0
    for (a_,b_) in cands:
        for wv in [(a_,b_),(-b_,a_)]:
            t2=nrm(wv)
            if t2==0 or t2<h1sq or t2*10**6>hU2: continue
            dd=det(wv,U)
            # t^2 |sin| = t^2 |det|/(t h_U) <1  <=> t^2 det^2 < hU2
            if t2*dd*dd<hU2:
                cnt+=1
                if not any(t2%l==0 for l in mods): ok=False; print("FAIL cand",wv)
    if trial<3: print("trial",trial,"close approximants tested",cnt)
print("ALL OK" if ok else "FAILURES")
