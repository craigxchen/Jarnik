# y_i=1/x_i: e1(y)=0, e2(y)e5(y)+e3(y)e4(y)=0, -e4(y) a nonzero square; then x_i = M/y_i
import itertools, sys
from math import gcd, isqrt
from functools import reduce
R=int(sys.argv[1]); out=[]
def esym(v):
    e=[1,0,0,0,0,0]
    for a in v:
        for k in range(5,0,-1): e[k]+=e[k-1]*a
    return e
for c in itertools.combinations(range(-R,R+1),4):
    y5=-sum(c)
    if y5 in c or 0 in c or y5==0: continue
    if y5<c[3]: continue   # canonical: y5 largest
    v=c+(y5,)
    e=esym(v)
    if e[2]*e[5]+e[3]*e[4]!=0: continue
    t=-e[4]
    if t<=0 or isqrt(t)**2!=t: continue
    M=reduce(lambda a,b:a*b//gcd(a,b),[abs(a) for a in v])
    x=[M//a for a in v]
    g=reduce(gcd,[abs(a) for a in x]); x=[a//g for a in x]
    ex=esym(x); s1=ex[1]
    b=[a*bb*(s1-a-bb) for a,bb in itertools.combinations(x,2)]
    ok = 0 not in b and len(set(b))==10
    out.append((v,x,ok))
    print('y',v,'x',x,'s=',ex[1:],'b distinct&nonzero',ok,flush=True)
print('total',len(out))
