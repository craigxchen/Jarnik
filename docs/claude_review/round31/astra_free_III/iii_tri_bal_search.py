from fractions import Fraction as F
from itertools import combinations, product
from math import gcd, lcm
import sympy as sp


def make_data(x,m):
    if x==0 or x==1 or not m or m*m==-x: return
    y=x*(1-x)/(m*m+x); z=1-x-y; u=m*y
    if x*y*z==0 or len({x,y,z})<3: return
    sets=[set()]
    for a,b,c in ((x,y,z),(x,z,y),(y,z,x)):
        sets.append({(F(0),F(1)),(u,a-1),(u,b-1),(u*(1+c)/c,c)})
    roots=set.union(*sets)
    if len(roots)!=7 or any(not im or (re,-im) in roots for re,im in roots): return
    D=lcm(*(v.denominator for r in roots for v in r))
    words=[]
    for i,j in combinations(range(4),2):
        loc=[(int(D*re),int(D*im)) for re,im in sets[i]-sets[j]]
        loc +=[(int(D*re),-int(D*im)) for re,im in sets[j]-sets[i]]
        if len(loc)!=4: return
        poly=[(1,0)]
        for re,im in loc:
            tmp=[(0,0)]*(len(poly)+1)
            for k,(a,b) in enumerate(poly):
                tmp[k]=(tmp[k][0]-a*re+b*im,tmp[k][1]-a*im-b*re)
                tmp[k+1]=(tmp[k+1][0]+a,tmp[k+1][1]+b)
            poly=tmp
        assert all(b==0 for a,b in poly[1:])
        if not poly[0][1]: return
        words.append((tuple(a for a,b in poly),abs(poly[0][1])))
    return D,words


def vp(n,p,B=999):
    if n==0:return B
    n=abs(n); k=0
    while k<B and n%p==0:
        n//=p;k+=1
    return k


def evaluate(poly,n):
    y=0
    for a in reversed(poly):
        y=y*n+a
    return y


def local(p,polys,heights,need_even=True):
    maxv=max(heights)
    base=list(range(p)) if p<=199 else None
    if base is None:
        t=sp.symbols('t'); base={0,1}
        for f,h in zip(polys,heights):
            if h:
                expr=sum(a*t**i for i,a in enumerate(f))
                _,factors=sp.factor_list(expr,modulus=p)
                for ff,mult in factors:
                    fp=sp.Poly(ff,t,modulus=p)
                    if fp.degree()==1:
                        a,b=map(int,fp.all_coeffs())
                        base.add((-b*pow(a,-1,p))%p)
        # include a simultaneous nonroot even if 0,1 happen to be roots
        for n in range(2,25):
            if all(evaluate(f,n)%p for f in polys):
                base.add(n);break
    def recurse(n,power,j):
        low=[]; can_descend=False
        for f,h in zip(polys,heights):
            value=evaluate(f,n)
            r=vp(value,p,min(h,j))
            if r==j and j<h:
                low.append((0,h-j));can_descend=True
            else:low.append((h-r,h-r))
        mn=max(a for a,b in low);mx=min(b for a,b in low)
        if p==2 and need_even: mn=max(1,mn)
        if mn>mx:return
        if not can_descend:return n,power
        # Try n itself (fix all p-adic digits), then branch if necessary.
        ts=[h-vp(evaluate(f,n),p,h) for f,h in zip(polys,heights)]
        if len(set(ts))==1 and (p!=2 or not need_even or ts[0]>0):
            return n,p**maxv
        if p>1000: return # optional search skip
        for k in range(p):
            res=recurse(n+k*power,power*p,j+1)
            if res:return res
    for n in base:
        res=recurse(n,p,1)
        if res:return res


CYCLES=[[0,1,4,5],[0,2,3,5],[1,2,3,4]]  # four cross edges, remaining two a matching


def search(x,m):
    data=make_data(x,m)
    if not data:return
    d,words=data
    facts=[sp.factorint(b) for f,b in words]
    for E in CYCLES:
        primes=sorted(set.union(*(set(facts[i]) for i in E)))
        if 2 not in primes:continue
        # Prioritize primes whose edge valuations differ.
        primes.sort(key=lambda p:(p!=2, len({facts[i].get(p,0) for i in E})==1, p))
        loc={}
        for p in primes:
            ans=local(int(p),[words[i][0] for i in E],[facts[i].get(p,0) for i in E])
            if ans is None:break
            loc[int(p)]=ans
        else:return d,words,E,loc


if __name__=='__main__':
    import sys
    ps=[(F(-3),F(3))]
    for b in range(1,4):
        for d in range(1,4):
            for a in range(-8,9):
                for c in range(1,7):
                    pair=(F(a,b),F(c,d))
                    if pair not in ps:ps.append(pair)
    for j,(x,m) in enumerate(ps):
        result=search(x,m)
        if result:
            D,words,E,loc=result
            print(x,m,D,E,[b for f,b in words],loc,flush=True);break
        if j in (0,5,20,100):print('searched',j,flush=True)
