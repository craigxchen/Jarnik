# find a primitive polynomial of degree N over F_p (own code); print coefficients low->high (monic)
import sys, itertools, random
def factor(n):
    f=[];d=2
    while d*d<=n:
        while n%d==0: f.append(d); n//=d
        d+=1
    if n>1: f.append(n)
    return sorted(set(f))
def pmulmod(a,b,m,p):
    # polys as lists low->high, m monic
    N=len(m)-1
    r=[0]*(2*N)
    for i,x in enumerate(a):
        if x:
            for j,y in enumerate(b):
                r[i+j]=(r[i+j]+x*y)%p
    for d in range(len(r)-1,N-1,-1):
        c=r[d]
        if c:
            for j in range(N+1):
                r[d-N+j]=(r[d-N+j]-c*m[j])%p
    return r[:N]
def ppow(e,m,p):
    N=len(m)-1
    res=[1]+[0]*(N-1); base=[0,1]+[0]*(N-2) if N>1 else [m[0]*-1%p]
    while e:
        if e&1: res=pmulmod(res,base,m,p)
        base=pmulmod(base,base,m,p); e>>=1
    return res
def is_primitive(m,p):
    N=len(m)-1; Q=p**N
    one=[1]+[0]*(N-1)
    if ppow(Q-1,m,p)!=one: return False
    for f in factor(Q-1):
        if ppow((Q-1)//f,m,p)==one: return False
    return True
p=int(sys.argv[1]); N=int(sys.argv[2])
random.seed(1)
while True:
    m=[random.randrange(p) for _ in range(N)]+[1]
    if m[0]==0: continue
    if is_primitive(m,p):
        print(' '.join(map(str,m))); break
