# Independent exact check of Prop 2.2 / Prop 6.1 (Hilbert 90 model), written from scratch.
import random
from fractions import Fraction as Fr
from math import gcd
random.seed(20261028)
PR=[5,13,17,29,37,41,53,61,73,89,97]
def gmul(a,b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def gconj(a): return (a[0],-a[1])
def gpow(a,k):
    r=(1,0)
    for _ in range(k): r=gmul(r,a)
    return r
def gnorm(a): return a[0]*a[0]+a[1]*a[1]
def gdivexact(a,b):  # a/b in Z[i] assumed exact
    n=gnorm(b); c=gmul(a,gconj(b))
    assert c[0]%n==0 and c[1]%n==0; return (c[0]//n,c[1]//n)
def gdivides(b,a):
    n=gnorm(b); c=gmul(a,gconj(b)); return c[0]%n==0 and c[1]%n==0
def gprime(p):
    for a in range(1,p):
        b2=p-a*a
        if b2<=0: break
        b=int(round(b2**0.5))
        if b*b==b2: return (a,b)
def val(a,pi):
    k=0
    while a!=(0,0) and gdivides(pi,a): a=gdivexact(a,pi); k+=1
    return k
UNITS=[(1,0),(0,1),(-1,0),(0,-1)]
def bracket(P,Q): return P[0]*Q[1]-P[1]*Q[0]   # Im(conj(P) Q)
stats=dict(trials=0,pairs=0,same=0,diff=0,fail=0,odd_half=0)
for trial in range(400):
    ps=random.sample(PR,random.randint(1,4))
    k={p:random.randint(1,3) for p in ps}
    pis={p:gprime(p) for p in ps}
    while True:
        A=[{p:random.randint(0,k[p]) for p in ps} for _ in range(5)]
        if all(min(a[p] for a in A)==0 and max(a[p] for a in A)==k[p] for p in ps): break
    us=[random.choice(UNITS) for _ in range(5)]
    Z=[]
    for j in range(5):
        z=us[j]
        for p in ps: z=gmul(z,gmul(gpow(pis[p],A[j][p]),gpow(gconj(pis[p]),k[p]-A[j][p])))
        Z.append(z)
    if len(set(Z))<5: continue
    N=gnorm(Z[0]); assert all(gnorm(z)==N for z in Z)
    stats['trials']+=1
    # x_j = Z_j/Z_1 ; P_j = primitive part of (Z_j*conj(Z_1) + N) [= N(1+x_j)], or i if x_j=-1
    Ps=[]
    for j in range(5):
        w=gmul(Z[j],gconj(Z[0])); P=(w[0]+N,w[1])
        if P==(0,0): P=(0,1)
        g=gcd(P[0],P[1]); P=(P[0]//g,P[1]//g)
        # check P/conj(P) = x_j  <=>  P*conj(Z_1)... : P*Z_1 conj?  x_j = Z_j/Z_1 ; P/Pbar = Z_j/Z_1 <=> P*Z_1 = Z_j*Pbar
        assert gmul(P,Z[0])==gmul(Z[j],gconj(P))
        Ps.append(P)
    eps=[{p:val(P,pis[p])-val(P,gconj(pis[p])) for p in ps} for P in Ps]
    e2=[val(P,(1,1)) for P in Ps]
    ok=True
    for j in range(5):
        nn=2**e2[j]
        for p in ps: nn*=p**abs(eps[j][p])
        ok&= (gnorm(Ps[j])==nn)
        for p in ps: ok&= (eps[j][p]-eps[0][p]==A[j][p]-A[0][p])
    H=1
    for p in ps: H*=p**(max(e[p] for e in eps)-min(e[p] for e in eps))
    ok&=(H==N)
    for a in range(5):
        for b in range(a+1,5):
            stats['pairs']+=1
            br=bracket(Ps[a],Ps[b]); C=1; ed=1
            for p in ps:
                ea,eb=eps[a][p],eps[b][p]
                c=max(0,min(ea,eb),min(-ea,-eb)); C*=p**c; ed*=p**abs(A[a][p]-A[b][p])
                ok&= (2*c-(abs(ea)+abs(eb))==-abs(ea-eb))
            ok&=(br%C==0) and br!=0
            # t^2 from the record's definition |z_a-z_b| = 2|t| R'' e^{-d/2}: t^2 = |z_a-z_b|^2 e^{d}/(4N)
            d=(Z[a][0]-Z[b][0],Z[a][1]-Z[b][1])
            t2=Fr(gnorm(d)*ed,4*N)
            t2b=Fr(br*br, C*C*2**(e2[a]+e2[b]))
            ok&=(t2==t2b)
            # unit class: e_a==e_b iff u_a/u_b in {+1,-1}
            ua=UNITS.index(us[a]); ub=UNITS.index(us[b]); sameclass=((ua-ub)%2==0)
            ok&= (sameclass==(e2[a]==e2[b]))
            if e2[a]==e2[b]:
                stats['same']+=1; ok&=(t2.denominator==1)
            else:
                stats['diff']+=1; ok&=(t2.denominator==2 and t2.numerator%2==1)
                stats['odd_half']+=1
    if not ok: stats['fail']+=1
print(stats)
