# Exact checks: (a) mnn31 points lie on the cone (D2-D1)(D3-D1)=(D4-D1)(D5-D1); (b) m-czech W points lie on none of the
# 15 S_5-images of the cone; (c) triple cross-product identity x^tau x y^tau = t^tau on both families (from the Cox data);
# (d) the twisted Galois involution u -> 1/u does not preserve the cone (symbolic).
import itertools, sympy as sp
from math import gcd
from sympy import factorint
def cone_forms(v):
    out=[]
    for base in range(5):
        rest=[i for i in range(5) if i!=base]
        for pairing in [((rest[0],rest[1]),(rest[2],rest[3])),((rest[0],rest[2]),(rest[1],rest[3])),((rest[0],rest[3]),(rest[1],rest[2]))]:
            (a,b),(c,d)=pairing
            out.append(((v[a]-v[base])*(v[b]-v[base])-(v[c]-v[base])*(v[d]-v[base]),base,pairing))
    return out
# (a) mnn31 Pell points
u,w=1843,1064
for k in range(6):
    D=[(u+2)*(u+3),(u+2)*(u+5),(u+3)*(u+5),3*(w+1)*(w+3),3*(w+2)*(w+3)]
    val=(D[1]-D[0])*(D[2]-D[0])-(D[3]-D[0])*(D[4]-D[0])
    assert val==0, val
    u,w=2*u+3*w+7,u+2*w+4
print("(a) cone identity exact at 6 mnn31 Pell points: OK")
# (b) m-czech W points: x on z^2=56x^2+65
pts=[]; x,z=1,11
# generate Pell solutions of z^2-56x^2=65 by brute small search then unit (15+2sqrt56)
sols=[(x,z) for x in range(1,200) for z in [sp.sqrt(56*x*x+65)] if z.is_integer]
print("small solutions:",sols[:5])
x0,z0=sols[0]
def nxt(x,z): return (15*x+2*z, 15*z+112*x)   # multiply by 15+2sqrt56 (norm 225-224=1)
cnt=0
for (x,z) in [(37,277)]:
    X,Z=x,z
    for k in range(5):
        X,Z=int(X),int(Z)
        S=(1+4+16+49)//1; m=1*2*4*7
        a,b,c=1,2,4; d=7
        F0=(X-a)*(X-b)*(X-c)*(X+d)
        e2=(a*a*b*b+a*a*c*c+a*a*d*d+b*b*c*c+b*b*d*d+c*c*d*d)
        F1=X**4-(a*a+b*b+c*c+d*d)*X*X//2+m+2*X*Z  # S/2 with S = sum of squares = 70
        F2=(X*X-b*b)*(X*X-c*c); F3=(X*X-a*a)*(X*X-c*c); F4=(X*X-a*a)*(X*X-b*b)
        n=(X*X-1)*(X*X-4)*(X*X-16)*(X*X-49)
        F=[F0,F1,F2,F3,F4]
        if 0 not in F and not all(n%f==0 for f in F):
            F[1]=X**4-70*X*X//2+m-2*X*Z
        if 0 not in F and all(n%f==0 for f in F) and len(set(F))==5:
            vals=[v for v,_,_ in cone_forms(F)]
            assert all(v!=0 for v in vals); cnt+=1
        X,Z=nxt(X,Z)
print("(b) W (m-czech) points checked:",cnt," none on any of the 15 cones")
# (c) triple cross product from Cox data: x^tau_a = prod_{B cap tau = {a}} n_B, y^tau_a = prod_{B cap tau = tau-{a}} n_B
def layers(u):
    g=0
    for t in u: g=gcd(g,t)
    u=[t//g for t in u]
    nB={}
    primes=set()
    for t in u: primes|=set(factorint(t).keys())
    for p in primes:
        v=[0]*5
        for i,t in enumerate(u):
            while t%p==0: t//=p; v[i]+=1
        for s in range(1,max(v)+1):
            B=frozenset(i for i in range(5) if v[i]>=s)
            nB[B]=nB.get(B,1)*p
    return u,nB
def triple_check(u):
    u,nB=layers(u)
    ok=True
    for tau in itertools.combinations(range(5),3):
        T=set(tau)
        X=[1,1,1];Y=[1,1,1]
        for B,val in nB.items():
            I=B&T
            if len(I)==1: X[tau.index(min(I))]*=val
            if len(I)==2: Y[tau.index(min(T-I))]*=val
        def tt(i,j): return (u[i]-u[j])//gcd(u[i],u[j])
        a,b,c=tau
        cross=(X[1]*Y[2]-X[2]*Y[1], X[2]*Y[0]-X[0]*Y[2], X[0]*Y[1]-X[1]*Y[0])
        tvec=(tt(b,c),tt(c,a),tt(a,b))
        ok&= (cross==tvec)
    return ok
u,w=1843,1064
D=[(u+2)*(u+3),(u+2)*(u+5),(u+3)*(u+5),3*(w+1)*(w+3),3*(w+2)*(w+3)]
print("(c) triple identity on mnn31 point:",triple_check(D))
import random
random.seed(5); okall=True
for _ in range(200):
    v=[random.randint(1,10**6) for _ in range(5)]
    if len(set(v))<5: continue
    okall&=triple_check(v)
print("(c) triple identity on 200 random 5-tuples:",okall)
# (d) Galois: u -> 1/u
uu=sp.symbols('u1:6')
Q=(uu[1]-uu[0])*(uu[2]-uu[0])-(uu[3]-uu[0])*(uu[4]-uu[0])
Qs=sp.together(Q.subs({uu[i]:1/uu[i] for i in range(5)}))
num=sp.factor(sp.numer(Qs))
print("(d) image of cone under u->1/u: numerator",num)
print("    divisible by Q?", sp.rem(sp.Poly(sp.numer(Qs),*uu),sp.Poly(Q,*uu)).is_zero)
