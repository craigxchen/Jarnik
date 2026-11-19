# Check, on random hyperplane-side families for random V (containing 1) and random T (small and large):
#  (a) F_i - F_j = P_J(1) P_J(Y+a)   (b) delta_ij := deg(F_i-F_j) - deg gcd(F_i,F_j) >= iota(P_J(1)) (odd-degree-prime part)
#  (c) the margin identity sum_{i<j} (D'/2 - deg gcd) <= s D'/4 for the reduced family (holds for ANY family with common norm,
#      cluster or not, as an identity on gcds: sum_{i<j} deg gcd >= D' s(s-2)/4)
#  (d) v_P(P_J(1)) >= |J cap K_P| unless one of the two forms vanishes on K_P = ker(ev_P|V).
import random, itertools, sys
from f2lib import *
def irreducibles(maxdeg):
    irr=[]
    for p in range(2,1<<(maxdeg+1)):
        if all(divmod2(p,q)[1] for q in irr if 2*deg(q)<=deg(p)): irr.append(p)
    return irr
IRR=irreducibles(10)
def vP(x,P):
    k=0
    while True:
        q,r=divmod2(x,P)
        if r: return k
        x=q;k+=1
def iota(x):  # degree of the part of x composed of odd-degree irreducibles (all factors here have degree <=10)
    s=0
    for P in IRR:
        if deg(P)%2==1:
            s+=deg(P)*vP(x,P)
    return s
def span(basis):
    S=[0]
    for b in basis: S+= [x^b for x in S]
    return S
def par(x): return bin(x).count('1')&1
if __name__!="__main__": raise SystemExit
random.seed(int(sys.argv[1])); trials=int(sys.argv[2])
nviol=0; cnt=0
for tr in range(trials):
    n=random.randint(3,5); e=random.randint(n,8)
    while True:
        basis=[1]+[random.randrange(2,1<<e) for _ in range(n-1)]
        V=span(basis)
        if len(set(V))==len(V): break
    # forms: coordinate vectors c in F_2^n with c_0=1 (l(1)=1); l(sum x_k b_k)=sum c_k x_k
    coords={}
    for x in range(1<<n):
        v=0
        for k in range(n):
            if x>>k&1: v^=basis[k]
        coords[v]=x
    forms=[c for c in range(1<<n) if c&1]
    s=random.randint(2,len(forms)); F=random.sample(forms,s); sides=[random.randint(0,1) for _ in F]
    T=random.choice([random.randrange(0,1<<random.randint(1,10)), (1<<random.randint(e,14))^random.randrange(1<<e)])
    Y=(T,1)
    S=[[v for v in V if par(coords[v]&c)==sd] for c,sd in zip(F,sides)]
    Fp=[prod4([add4(Y,(v,0)) for v in Si]) for Si in S]
    G=Fp[0]
    for f in Fp[1:]: G=gcd4(G,f)
    D0=deg4(Fp[0]); Dp=D0-deg4(G)
    sumg=0
    for i,j in itertools.combinations(range(s),2):
        J=[v for v in V if par(coords[v]&F[i])==0 and par(coords[v]&F[j])==0]
        a=[v for v in S[i] if v in set(S[j])][0]
        PJ1=1
        for v in J: PJ1=mul(PJ1,1^v)
        PJYa=prod4([add4(Y,(a^v,0)) for v in J])
        df=add4(Fp[i],Fp[j]); assert df==m4((PJ1,0),PJYa), "identity (a)"
        g=gcd4(Fp[i],Fp[j]); delta=deg4(df)-deg4(g)
        sumg+=deg4(g)-deg4(G)
        io=iota(PJ1); cnt+=1
        if delta<io: nviol+=1; print("VIOLATION b",basis,T,F[i],F[j],delta,io)
        # (d)
        for P in IRR:
            if deg(P)%2==0 or deg(P)>e: continue
            K=[v for v in V if divmod2(v,P)[1]==0]
            spec=lambda c: all(par(coords[k]&c)==0 for k in K)
            if not spec(F[i]) and not spec(F[j]):
                JK=[v for v in J if v in set(K)]
                assert vP(PJ1,P)>=len(JK), "(d)"
    assert 4*sumg>=Dp*s*(s-2), "(c) margin identity"
print("pairs checked",cnt,"violations of (b)",nviol)
