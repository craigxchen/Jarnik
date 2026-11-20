import random, itertools
from math import prod
from sympy import primerange
from circle_sim import two_squares, gmul
rng=random.Random(5)
pool=[p for p in primerange(100,2000) if p%4==1]
def conj(g): return (g[0],-g[1])
def gp(gs): 
    z=(1,0)
    for g in gs: z=gmul(z,g)
    return z
cnt=0
for trial in range(300):
    pr=rng.sample(pool,10); gens=[two_squares(p) for p in pr]
    eta=rng.choice([(1,0),(0,1),(-1,0),(0,-1)])
    ch=[[rng.randrange(2) for _ in pr] for _ in range(3)]
    if len({tuple(c) for c in ch})<3: continue
    z=[gmul(eta,gp([g if c==0 else conj(g) for g,c in zip(gens,cc)])) for cc in ch]
    A={'12':[], '23':[], '13':[]}; G=[]
    for k,g in enumerate(gens):
        b=[ch[i][k] for i in range(3)]
        f=g if b[0]==0 else conj(g)
        if b[0]==b[1]==b[2]: G.append(f)
        elif b[0]==b[1]: A['12'].append(f)
        elif b[1]==b[2]: A['23'].append(g if b[1]==0 else conj(g))
        else: A['13'].append(f)
    def tval(i,j,key):
        H=gp(G+A[key]); d=(z[j][0]-z[i][0],z[j][1]-z[i][1])
        den=gmul(gmul((0,2),eta),H); nd=den[0]**2+den[1]**2
        num=gmul(d,conj(den)); assert num[0]%nd==0 and num[1]%nd==0
        t=(num[0]//nd,num[1]//nd); assert t[1]==0 and t[0]!=0; return t[0]
    t12,t23,t13=tval(0,1,'12'),tval(1,2,'23'),tval(0,2,'13')
    a12,a23,a13=gp(A['12']),gp(A['23']),gp(A['13'])
    lhs=(t12*a12[0]+t23*a23[0],t12*a12[1]+t23*a23[1]); rhs=(t13*a13[0],t13*a13[1])
    assert lhs==rhs
    # norm of the forced block determines its rational primes
    for key in ('13','23'):
        assert gp(A[key])[0]**2+gp(A[key])[1]**2==prod(g[0]**2+g[1]**2 for g in A[key])
    cnt+=1
print('Gaussian triangle identity verified on',cnt,'random triples')
