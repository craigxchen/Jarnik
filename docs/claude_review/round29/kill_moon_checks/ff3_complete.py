"""Complete F_3(t) scan for D<=DMAX: ALL squarefree split N (prime pairs of every degree <= DMAX), via python-flint
factorisation of F_3-irreducibles of degree 2d over F_9. Reuses arithmetic of ff3_scan.py."""
import sys, itertools, flint
sys.argv=[sys.argv[0]]+sys.argv[1:]
DMAX=int(sys.argv[1])
import importlib.util
src=open('ff3_scan.py').read().split('DMAX=int')[0]   # field + poly code only
exec(src)
m=flint.fmpz_mod_poly_ctx(3)([1,0,1]); F=flint.fq_default_ctx(modulus=m,var='i'); R9=flint.fq_default_poly_ctx(F)
def conv(f):
    out=[]
    for c in f.coeffs():
        L=c.to_list()+[0,0]; out.append(idx[(int(L[0])%3,int(L[1])%3)])
    return tuple(out)
pairs={}
for d in range(1,DMAX+1):
    L=[]
    for co in itertools.product(range(3),repeat=2*d):
        P=flint.fmpz_mod_poly_ctx(3)(list(co)+[1])
        if not P.is_irreducible(): continue
        fac=R9([int(x) for x in list(co)+[1]]).factor()[1]
        assert len(fac)==2
        L.append((conv(fac[0][0]),conv(fac[1][0])))
    pairs[d]=L
print({d:len(v) for d,v in pairs.items()})
allp=[(d,p) for d in sorted(pairs) for p in pairs[d]]
from collections import defaultdict
best=defaultdict(int); nN=defaultdict(int); arg=defaultdict(lambda:None)
def evaluate(ch,D):
    nN[D]+=1
    full=[]
    for eps in itertools.product((0,1),repeat=len(ch)):
        z=(ONE,)
        for (d,(f,g)),e in zip(ch,eps): z=pmul(z,f if e==0 else g)
        full.append(z)
    for c in (-1,0,1):
        s=int((D/2+c)//1)
        groups=defaultdict(int)
        for z in full: groups[z[s+1:]]+=1
        mx=max(groups.values())
        if mx>best[(D,c)]: best[(D,c)]=mx; arg[(D,c)]=[d for d,_ in ch]
def rec(start,D,chosen):
    for k in range(start,len(allp)):
        d,p=allp[k]
        if D+d>DMAX: break
        ch=chosen+[(d,p)]; evaluate(ch,D+d); rec(k+1,D+d,ch)
rec(0,0,[])
print("N counts by D:",dict(nN))
for D in range(1,DMAX+1): print("D=",D," max cluster c=-1,0,1:",[best[(D,c)] for c in (-1,0,1)], "attained with pair degrees", [arg[(D,c)] for c in (-1,0,1)])
