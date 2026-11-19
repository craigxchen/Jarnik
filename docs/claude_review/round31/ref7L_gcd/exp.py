import sys, itertools, random, math
sys.path.insert(0,'.')
from gcdlib import *
random.seed(1)
def rnd(): return random.randrange(1,P)
box=[w for w in itertools.product(range(-2,3),repeat=3) if any(w)]
def bestpair(n):
    best=(0,None,None)
    vals={w:sum(a*b for a,b in zip(w,n)) for w in box}
    for i,w in enumerate(box):
        for w2 in box[i+1:]:
            # independent
            c=[w[1]*w2[2]-w[2]*w2[1],w[2]*w2[0]-w[0]*w2[2],w[0]*w2[1]-w[1]*w2[0]]
            if not any(c): continue
            a,b=vals[w],vals[w2]
            g=math.gcd(a,b)
            if a==0 and b==0: g=10**9
            if g>best[0]: best=(g,w,w2)
    return best
def mono(w,c=1): return {w:c}
def add(*Fs):
    r={}
    for F in Fs:
        for u,c in F.items(): r[u]=(r.get(u,0)+c)%P
    return {u:c for u,c in r.items() if c}
def mul(F,G):
    r={}
    for u,a in F.items():
        for v,b in G.items():
            k=tuple(x+y for x,y in zip(u,v)); r[k]=(r.get(k,0)+a*b)%P
    return {u:c for u,c in r.items() if c}
def binom(w): return add({w:1},{(0,0,0):P-1})
def rpoly(): # 2 random monomials in [0,1]^3, random coefs
    S=random.sample(list(itertools.product(range(2),repeat=3)),2)
    return {u:rnd() for u in S}
seqs={
 'generic (l, l*sqrt2, l*sqrt3)': lambda l:(l,round(l*2**.5),round(l*3**.5)),
 'adversarial (m, m*p, m*p+1), p=floor(log2 l)+1': None,
}
def adv(L):
    p=max(2,int(math.log2(L)))
    m=L//p
    return (m, m*p-1 if False else m*p, m*p - m//2 + 1), p
print("T2/T3: worst binomial-pair gcd over w,w' in [-2,2]^3, and planted random pair (mod P gcd)")
for name in ['gen','adv']:
    for L in [64,128,256,512,1024,2048]:
        if name=='gen':
            l=round(L/3**.5); n=(l,round(l*2**.5),round(l*3**.5))
        else:
            n,p=adv(L)
        g,w,w2=bestpair(n)
        Bw,Bw2=binom(w),binom(w2)
        F=add(mul(rpoly(),Bw),mul(rpoly(),Bw2)); G=add(mul(rpoly(),Bw),mul(rpoly(),Bw2))
        d=common([F,G],n)
        Lr=max(abs(x) for x in n)
        print(name,n,'L=',Lr,'best binom gcd',g,'w',w,w2,'planted gcddeg',d,'ratio',round(d/Lr,4))
