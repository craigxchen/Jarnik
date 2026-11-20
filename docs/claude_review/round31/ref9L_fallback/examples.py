import itertools, random, math
from sympy import primerange, factorint
from sympy.ntheory.modular import crt
Ts=[T for r in range(3,7) for T in itertools.combinations(range(7),r)]
pairs=list(itertools.combinations(range(7),2))
def cof(y,e):
    c={}
    for i,j in pairs:
        pr=1
        for T in Ts:
            if i in T and j in T: pr*=e[T]
        assert (y[i]-y[j])%pr==0
        c[(i,j)]=(y[i]-y[j])//pr
    return c
def bad(c,exempt):
    # primes p not in exempt dividing >=2 cofactors
    cnt={}
    for v in c.values():
        for p in factorint(abs(v), limit=10**6):
            if p not in exempt: cnt[p]=cnt.get(p,0)+1
    return {p:k for p,k in cnt.items() if k>=2}
# (1) census-style: 98 distinct primes ~ A for |T|>=3, random lattice vector by CRT
random.seed(2)
A=1000
pr=list(primerange(A,3*A))[:98]; e={T:pr[k] for k,T in enumerate(Ts)}
for trial in range(5):
    y=[]
    res={T:random.randrange(e[T]) for T in Ts}
    for i in range(7):
        mods=[e[T] for T in Ts if i in T]; rs=[res[T] for T in Ts if i in T]
        r,M=crt(mods,rs); y.append(int(r)+random.randrange(5)*int(M))
    c=cof(y,e)
    ev=sum(1 for v in c.values() if v%2==0); t3=sum(1 for v in c.values() if v%3==0)
    print("census-style vector: even cofactors",ev,"multiples of 3",t3)
# (2) explicit integer input for the corrected local conditions (no balance): e_{0,1,2}=5, other e_T=1
e1={T:1 for T in Ts}; e1[(0,1,2)]=5
found=0
for trial in range(200000):
    y=[0,5,10]+[random.randrange(-400,400) for _ in range(4)]
    if len(set(y))<7: continue
    g=0
    for v in y: g=math.gcd(g,v)
    if g!=1: continue
    c=cof(y,e1)
    if not bad(c,{2,3}):
        print("corrected-input example y=",y,"cofactors",[c[p] for p in pairs]); found+=1
        if found>=2: break
