"""Independent reproduction of the fourteen-prime computations of III.2 / III.7 (answer.md lines 730, 898-902).
Own solver: exact kernel basis of (A,B,C), Lagrange reduction, sup-norm enumeration."""
import sys
from itertools import combinations
from math import gcd, isqrt, prod

def isprime(n):
    if n<2: return False
    for p in (2,3,5,7,11,13,17,19,23,29,31,37):
        if n%p==0: return n==p
    d=n-1; s=0
    while d%2==0: d//=2; s+=1
    for a in (2,3,5,7,11,13,17,19,23,29,31,37):
        x=pow(a,d,n)
        if x in (1,n-1): continue
        for _ in range(s-1):
            x=x*x%n
            if x==n-1: break
        else: return False
    return True

def egcd(a,b):
    if b==0: return (a,1,0)
    g,x,y=egcd(b,a%b); return (g,y,x-(a//b)*y)

def dot(u,v): return sum(a*b for a,b in zip(u,v))

def solutions(A,B,C,U,first_only=False):
    """All primitive (up to sign) integer (x,y,z), all nonzero, |.|<=U, with A x + B y + C z = 0."""
    g,a,b=egcd(A,B); assert g==1
    v=(B,-A,0); w=(-C*a,-C*b,1)
    assert dot(v,(A,B,C))==0 and dot(w,(A,B,C))==0
    while True:  # Lagrange
        if dot(v,v)>dot(w,w): v,w=w,v
        k=(2*dot(v,w)+dot(v,v))//(2*dot(v,v))
        if k==0: break
        w=tuple(y-k*x for x,y in zip(v,w))
    if dot(v,v)>3*U*U: return []
    vv,ww,vw=dot(v,v),dot(w,w),dot(v,w); D=vv*ww-vw*vw   # D = |w*|^2 vv
    # |beta|*|w*| <= |vec| <= sqrt(3) U  ->  beta^2 * D <= 3U^2 vv
    bmax=isqrt(3*U*U*vv//D)+1
    out=[]
    for beta in range(0,bmax+1):
        lo,hi=-10**40,10**40
        for k in range(3):
            c=v[k]; e=beta*w[k]
            # |alpha*c+e|<=U
            if c==0:
                if abs(e)>U: lo,hi=1,0
                continue
            if c>0: l=-((U+e)//c); h=(U-e)//c   # ceil((-U-e)/c), floor((U-e)/c)
            else: l=-((U-e)//(-c)); h=(U+e)//(-c)
            lo=max(lo,l); hi=min(hi,h)
        for alpha in range(lo,hi+1):
            if beta==0 and alpha<=0: continue
            vec=tuple(alpha*x+beta*y for x,y in zip(v,w))
            if all(vec) and max(map(abs,vec))<=U and gcd(*vec)==1:
                assert dot(vec,(A,B,C))==0
                out.append(vec)
                if first_only: return out
    return out

def partitions(n,k):
    """unordered triples of disjoint k-subsets of range(n)"""
    subs=list(combinations(range(n),k)); idx={s:i for i,s in enumerate(subs)}
    for i,a in enumerate(subs):
        for j in range(i+1,len(subs)):
            b=subs[j]
            if set(a)&set(b): continue
            for l in range(j+1,len(subs)):
                c=subs[l]
                if set(c)&(set(a)|set(b)): continue
                yield a,b,c

def run(ps,k,label):
    N=prod(ps); s=sorted(ps)
    den=16*isqrt(N); U=-(-prod(s[-8:])//den)
    tot=0; passing=[]; indiv=prim=0; e=isqrt(U)
    for a,b,c in partitions(14,k):
        tot+=1
        A,B,C=(prod(ps[i] for i in t) for t in (a,b,c))
        sol=solutions(A,B,C,U)
        if sol:
            passing.append(((a,b,c),sol))
    # factorability: coefficient = product of two integers <= isqrt(U) (resp. two even after scaling by 4)
    def splits(v,bd):
        v=abs(v)
        return any(v%d==0 and v//d<=bd for d in range(1,min(bd,isqrt(v))+1))
    for _,sol in passing:
        # all solutions are multiples of the primitive one when det > 3U^2 (checked)
        assert len(sol)==1
        indiv+=all(splits(x,e) for x in sol[0]); prim+=all(splits(x,e//2) for x in sol[0])
    print(f"{label}: k={k} U={U} triples={tot} passing={len(passing)} indiv-factorable={indiv} even-factorable={prim}",flush=True)
    return U,tot,passing

if __name__=='__main__':
    which=sys.argv[1]
    if which=='first':
        ps=[]; c=100001
        while len(ps)<14:
            if c%4==1 and isprime(c): ps.append(c)
            c+=1
        print('first parent primes',ps[0],'...',ps[-1])
        run(ps,2,'first parent'); run(ps,3,'first parent')
    else:
        PL=[1000000873621,1000000976641,1000001041141,1000000874941,1000000973341,1000001043121,
            1000000876981,1000000968661,1000001045761,1000000963813,1000000963837,1000000963877,1000000963937,1000000964021]
        assert all(p%4==1 and isprime(p) for p in PL) and len(set(PL))==14
        ps=sorted(PL)
        U,tot,passing=run(ps,3,'planted parent')
        for (a,b,c),sol in passing:
            print(' passing partition (primes):',[[ps[i] for i in t] for t in (a,b,c)],'solution',sol)
