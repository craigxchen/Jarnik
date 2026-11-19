import numpy as np, itertools, math
from fractions import Fraction as Fr
rng=np.random.default_rng(7)
# 1. GMD <= (2/sqrt3) sd on random discrete laws (and equality check on uniform)
worst=0
for t in range(20000):
    k=rng.integers(2,12); x=rng.integers(-20,20,size=k).astype(float); p=rng.random(k)**3; p/=p.sum()
    gmd=sum(p[i]*p[j]*abs(x[i]-x[j]) for i in range(k) for j in range(k))
    mu=(p*x).sum(); sd=math.sqrt((p*(x-mu)**2).sum())
    if sd>0: worst=max(worst,gmd/sd)
print("max GMD/sd over random laws:",worst,"bound 2/sqrt3=",2/math.sqrt(3))
# 2. centered square mean over nonzero residues = q(q+1)/12
for q in (3,5,7,11,13):
    c=[(x if x<=q//2 else x-q) for x in range(1,q)]
    print(q, Fr(sum(v*v for v in c),q-1), Fr(q*(q+1),12))
# 3. numeric margins
eta=Fr(1,200); r=8*(1+eta)/(15*(1-eta)); print("8(1+eta)/(15(1-eta))=",float(r),"<13/24=",float(Fr(13,24)), r<Fr(13,24), "G>=",float(1-2*r),">= -1/12:",1-2*r>=-Fr(1,12))
for q in (3,5,7,11,101):
    print(q, Fr(4,9)*Fr(5*q*(q+1),12) <= Fr(20,81)*q*q, float(Fr(5*q*(q+1),27)/q**2))
print("k bound ~", 1/(1-2*math.sqrt(20/81)))
# 4. spectral lemma: G PSD with nonpositive off-diagonal & unit diagonal => sum_{j!=i} G_ij^2 <= 1
mx=0
for t in range(3000):
    n=rng.integers(3,30); A=rng.random((n,n))**rng.integers(1,6); A=(A+A.T)/2; np.fill_diagonal(A,0)
    lam=np.linalg.eigvalsh(A).max(); A/=lam*(1+1e-12)  # now lambda_max(A)=1 => G=I-A PSD
    G=np.eye(n)-A
    assert np.linalg.eigvalsh(G).min()>-1e-9
    mx=max(mx,(A**2).sum(axis=1).max())
print("max row sum of G_ij^2 (j!=i):",mx)
# 5. variance identity under layer+flip: Var Z = c^T G c / 4, on a random integer profile
m,P=5,6; e=rng.integers(1,4,size=P); logp=np.log(rng.choice([5,13,17,29,37,41],size=P,replace=False))
a=np.array([[rng.integers(0,ei+1) for ei in e] for _ in range(m)])
c=np.array([2,-1,1,-3,1])
layers=[(p,h) for p in range(P) for h in range(1,e[p]+1)]; W=sum(e*logp)
wts=np.array([logp[p] for p,h in layers])/W
B=np.array([[1 if a[i,p]>=h else 0 for (p,h) in layers] for i in range(m)])
f=2*B-1; G=(f*wts)@f.T
Z=c@B; Zf=c.sum()-Z
mean=0.5*(wts@Z+wts@Zf); var=0.5*(wts@(Z-mean)**2+wts@(Zf-mean)**2)
print("Var Z (layer+flip) =",var," c^TGc/4 =",c@G@c/4, " mean =",mean, "sum c/2=",c.sum()/2)
print("G_ij = 1-2w(a_i-a_j)/W check:", np.allclose(G, 1-2*np.array([[np.sum(np.abs(a[i]-a[j])*logp) for j in range(m)] for i in range(m)])/W))
