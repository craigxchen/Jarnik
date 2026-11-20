import numpy as np, itertools, random
exec(open('paley_half.py').read().split('for m in')[0])
# (5) exponents of the Paley character on circles
for m in [11,19,43]:
    s=S(m); b=(s+1)//2
    ok=True
    for v,w in itertools.combinations(range(m),2):
        a=s[:,v]-s[:,w]
        e=a@b
        exp=np.zeros(m,dtype=int); exp[v]=(m+1)//2; exp[w]=-(m+1)//2
        ok&= (e==exp).all() and a.sum()==0
    print("m",m,"a=S_v-S_w: sum a=0 and a.b_h = +-(m+1)/2 on h=v,w, 0 elsewhere:",ok, " ||a||_1 max", max(np.abs(s[:,v]-s[:,w]).sum() for v,w in itertools.combinations(range(m),2)))
# (4) Gram entries on circles: z_i conj(z_j) = h * unit * X^2, X of norm Q=q/h
def gmul(a,b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def conj(a): return (a[0],-a[1])
from sympy import primerange
def gauss_prime(p):
    for x in range(1,int(p**.5)+1):
        y2=p-x*x; y=int(round(y2**.5))
        if y*y==y2: return (x,y)
pr=[p for p in primerange(5,400) if p%4==1][:8]; q=1
for p in pr: q*=p
pis=[gauss_prime(p) for p in pr]
random.seed(3); bad=0; tested=0
for _ in range(300):
    sig=[[random.randint(0,1) for _ in pr] for _ in range(3)]
    zs=[]
    for sg in sig:
        z=(1,0)
        for bit,pi in zip(sg,pis): z=gmul(z, pi if bit else conj(pi))
        zs.append(z)
    for i,j in itertools.combinations(range(3),2):
        G=gmul(zs[i],conj(zs[j]))
        h=1;Q=1
        for k,p in enumerate(pr):
            if sig[i][k]==sig[j][k]: h*=p
            else: Q*=p
        X=(1,0)
        for k,pi in enumerate(pis):
            if sig[i][k]!=sig[j][k]: X=gmul(X, pi if sig[i][k] else conj(pi))
        X2=gmul(X,X); tested+=1
        if not any(G==(h*u[0]*X2[0]-h*u[1]*X2[1], h*u[0]*X2[1]+h*u[1]*X2[0]) for u in [(1,0),(-1,0),(0,1),(0,-1)]): bad+=1
    # triangle Gram determinant (real parts) times q^3 vanishes
    R={}
    for i,j in itertools.combinations(range(3),2): R[i,j]=gmul(zs[i],conj(zs[j]))[0]
    d=q**3+2*R[0,1]*R[1,2]*R[0,2]-q*(R[0,1]**2+R[1,2]**2+R[0,2]**2)
    assert d==0
print("Gram entries z_i conj z_j = h*unit*X_ij^2 failures:",bad,"of",tested,"; q^3 det(Re-Gram)=0 on all 300 genuine triples")
