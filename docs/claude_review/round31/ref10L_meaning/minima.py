# Check: for genuine twin configurations, the integer vectors z in M_S (z mod p constant on each class)
# with spread <= L span, together with 1, a space of dimension <= 3 (k <= 2).
# Upper bound via LLL Gram-Schmidt: short vectors lie in span(b_1..b_T), T = max{t: ||b*_t|| <= R}.
import random, math
from fractions import Fraction
from sympy import primerange
from sympy.matrices.normalforms import hermite_normal_form
from sympy import Matrix

def lll(B, delta=Fraction(99,100)):
    B=[list(map(int,r)) for r in B]; n=len(B)
    def gs(B):
        Bs=[];mu=[[Fraction(0)]*n for _ in range(n)]
        for i in range(n):
            v=[Fraction(x) for x in B[i]]
            for j in range(i):
                mu[i][j]=sum(Fraction(a)*b for a,b in zip(B[i],Bs[j]))/sum(b*b for b in Bs[j])
                v=[a-mu[i][j]*b for a,b in zip(v,Bs[j])]
            Bs.append(v)
        return Bs,mu
    Bs,mu=gs(B);k=1
    while k<n:
        for j in range(k-1,-1,-1):
            r=round(mu[k][j])
            if r:
                B[k]=[a-r*b for a,b in zip(B[k],B[j])]
                Bs,mu=gs(B)
        nk=sum(x*x for x in Bs[k]); nk1=sum(x*x for x in Bs[k-1])
        if nk>=(delta-mu[k][k-1]**2)*nk1: k+=1
        else:
            B[k],B[k-1]=B[k-1],B[k]; Bs,mu=gs(B); k=max(k-1,1)
    return B,[math.sqrt(float(sum(x*x for x in v))) for v in Bs]

def basis_Nprime(y,primes):
    B=len(y); q=math.prod(primes); gens=[]
    for i in range(B-1):
        e=[0]*(B-1); e[i]=q; gens.append(e)
    for p in primes:
        c=q//p; u=pow(c,-1,p)
        gens.append([(c*u)*((y[i+1]-y[0])%p!=0) for i in range(B-1)])
    H=hermite_normal_form(Matrix(gens).T).T  # rows basis? normalise
    rows=[list(H.row(i)) for i in range(H.rows) if any(H.row(i))]
    return rows,q

random.seed(1); CMAX=6.0
found=0; worst=0
for trial in range(200000):
    B=random.choice([4,5,6]); X=random.randint(30,3000)
    y=sorted(random.sample(range(X+1),B))
    if y[0]!=0: y=[v-y[0] for v in y]
    spread=y[-1]
    pr=[p for p in primerange(3,spread+1) if len({v%p for v in y})==2]
    if not pr: continue
    q=math.prod(pr)
    C=spread/math.sqrt(q)
    if C>CMAX: continue
    # need distinct signatures? not needed for the lattice statement
    rows,q=basis_Nprime(y,pr)
    if len(rows)!=B-1: print("rank?",len(rows)); continue
    red,gsn=lll(rows)
    R=math.sqrt(B-1)*spread
    T=max([t+1 for t in range(B-1) if gsn[t]<=R],default=0)
    # n! V_n bound prediction
    found+=1; worst=max(worst,T)
    print(f"B={B} y={y} s={len(pr)} logq={math.log(q):.1f} C={C:.3f} GS={[round(g,1) for g in gsn]} R={R:.1f} k_upper={T}")
    if found>=40: break
print("cases",found,"max k upper bound",worst)
