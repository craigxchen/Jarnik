# Referee check: (D4) parity rule, (D8) cut identity, (M1) product identity, level capacity,
# deficit lemma, Theorem A chain, layer lemma (Theorem C), literal vs corrected Theorem C,
# 6 | L_0 for M>=5.  All exact except the final log comparisons (done with a margin).
import random, math, sys
from fractions import Fraction
from math import gcd
from rlib import *

random.seed(int(sys.argv[1]) if len(sys.argv)>1 else 1)
SPLIT=[5,13,17,29,37,41,53,61,73,89,97,101,109,113]
def lcm(a,b): return a*b//gcd(a,b)

stats=dict(clusters=0,D4=0,D8=0,M1=0,cap=0,deficit=0,thmA=0,layer=0,lit_fail=0,six=0,six_tests=0,thmA_tight=0)
worstA=0
for trial in range(4000):
    k=random.randint(1,4)
    ps=random.sample(SPLIT[:9],k)
    N0=1
    for p in ps: N0*=p**random.randint(1,4)
    if N0>3*10**9: continue
    # extra non-primitive factor sometimes
    if random.random()<0.3: N0*=random.choice([2,4,9,5,25])
    pts=points_on_circle(N0)
    if len(pts)<4: continue
    pts.sort(key=angle)
    M=random.randint(3,min(9,len(pts)))
    # consecutive window or random subset inside an open half circle
    if random.random()<0.6:
        st=random.randrange(len(pts)); sub=[pts[(st+t)%len(pts)] for t in range(M)]
    else:
        sub=random.sample(pts,M)
    # require all in an open half-plane: angular span < pi
    angs=sorted(angle(z) for z in sub)
    gaps=[(angs[(t+1)%M]-angs[t])%(2*math.pi) for t in range(M)]
    if max(gaps)<=math.pi+1e-12: continue
    P,g=primitive_cluster(sub)
    N=gnorm(P[0]); assert all(gnorm(z)==N for z in P)
    stats['clusters']+=1
    pairs=list(itertools.combinations(range(M),2))
    s={};n={};eps={}
    L0=1
    for (i,j) in pairs:
        c=cot_half(P[i],P[j],N); s[i,j]=c.denominator; L0=lcm(L0,s[i,j])
        n[i,j]=N//gnorm(ggcd(P[i],P[j]))
        eps[i,j]=2 if (P[i][0]-P[j][0])%2 else 1
        d2=gnorm((P[i][0]-P[j][0],P[i][1]-P[j][1]))
        if d2*eps[i,j]*n[i,j]==4*N*s[i,j]**2: stats['D4']+=1
        else: print('D4 FAIL',N,P[i],P[j]); 
    fN=factor(N)
    assert all(p%4==1 for p in fN), fN
    levels={}; Dp4={}
    for p,e in fN.items():
        pi=gauss_prime_above(p)
        a=[vpi(z,pi) for z in P]
        assert min(a)==0 and max(a)==e
        levels[p]=a
        Dp4[p]=sum((2*sum(1 for x in a if x<tau)-M)**2 for tau in range(1,e+1))  # = 4 D_p
    # (D8): prod n^4 * prod p^{4D_p} = N^{M^2}
    lhs=1
    for e_ in pairs: lhs*=n[e_]**4
    for p in fN: lhs*=p**Dp4[p]
    if lhs==N**(M*M): stats['D8']+=1
    else: print('D8 FAIL')
    # (M1) 4th power form
    Pp=len(pairs)
    prodd2=1; prode=1; prods=1
    for (i,j) in pairs:
        prodd2*=gnorm((P[i][0]-P[j][0],P[i][1]-P[j][1])); prode*=eps[i,j]; prods*=s[i,j]
    lhs=prodd2**4*prode**4*N**(M*M); rhs=4**(4*Pp)*N**(4*Pp)*prods**8
    for p in fN: rhs*=p**Dp4[p]
    if lhs==rhs: stats['M1']+=1
    else: print('M1 FAIL')
    # capacity and deficit, Theorem A
    W=math.log(N); C=math.sqrt(max(gnorm((P[i][0]-P[j][0],P[i][1]-P[j][1])) for (i,j) in pairs))/N**0.25
    lhsA=0.0; sumD=0.0
    ok=True
    for p,e in fN.items():
        v=0; t=L0
        while t%p==0: t//=p; v+=1
        c=(p-1)*p**v
        cnt=[levels[p].count(h) for h in range(e+1)]
        if max(cnt)>c: ok=False; print('CAP FAIL',p,cnt,c)
        D=Dp4[p]/4
        if M>=4*c:
            if D+1e-9>=(M-4*c)**3/(12*c): stats['deficit']+=1
            else: print('DEFICIT FAIL',M,c,cnt)
            lhsA+=math.log(p)*(M-4*c)**3/(12*c)
        sumD+=D*math.log(p)
    # also split primes not dividing N: must have M <= c_p
    for p in SPLIT:
        if p in fN: continue
        v=0; t=L0
        while t%p==0: t//=p; v+=1
        if M>(p-1)*p**v: ok=False; print('CAP FAIL (p not | N)',p,M,L0)
    if ok: stats['cap']+=1
    rhsA=M*W/4+Pp*math.log(C*C/2)
    if lhsA<=sumD+1e-9 and sumD<=rhsA+1e-9: stats['thmA']+=1
    else: print('THM A FAIL',lhsA,sumD,rhsA)
    worstA=max(worstA, sumD/rhsA if rhsA>0 else 0)
    # layer lemma
    prodz=(1,0)
    for z in P: prodz=gmul(prodz,z)
    m=1; normPhi=1; Nbal=1; Delta=1
    for p,e in fN.items():
        for tau in range(1,e+1):
            h=sum(1 for x in levels[p] if x>=tau)
            m*=p**min(h,M-h); normPhi*=p**abs(M-2*h)
            if 2*h==M: Nbal*=p
            Delta*=p**((2*h-M)**2//4) if M%2==0 else 1
    Phi=(prodz[0]//m, prodz[1]//m) if prodz[0]%m==0 and prodz[1]%m==0 else None
    if Phi is not None and gnorm(Phi)==normPhi: stats['layer']+=1
    else: print('LAYER FAIL')
    if M%2==0:
        # literal statement: prod z = unit * Nbal^{M/2} * Phi' with Norm(Phi') <= Delta^2 ?
        q=gdivmod_exact(prodz,(Nbal**(M//2),0))
        if q is None or gnorm(q)>Delta**2: stats['lit_fail']+=1
        else:
            assert normPhi<=Delta**2
    if M>=5:
        stats['six_tests']+=1
        if L0%6==0: stats['six']+=1
        else: print('6|L0 FAIL',N,P)
print(stats); print('max sumD/rhsA',worstA)
