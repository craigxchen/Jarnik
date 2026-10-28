# exact Gaussian-integer helpers; Gaussian integers as (a,b) tuples
from math import gcd, log
from itertools import combinations
def mul(x,y): return (x[0]*y[0]-x[1]*y[1], x[0]*y[1]+x[1]*y[0])
def conj(x): return (x[0],-x[1])
def norm(x): return x[0]*x[0]+x[1]*x[1]
def pw(x,k):
    r=(1,0)
    for _ in range(k): r=mul(r,x)
    return r
def rdiv(a,b):  # nearest-integer quotient a/b (exact rounding)
    n=norm(b); t=mul(a,conj(b))
    def rnd(p,q): return (2*p+q)//(2*q)
    return (rnd(t[0],n), rnd(t[1],n))
def sub(a,b): return (a[0]-b[0],a[1]-b[1])
def ggcd(a,b):
    while b!=(0,0):
        q=rdiv(a,b); a,b=b,sub(a,mul(q,b))
    return a
def exdiv(a,b):  # exact division
    n=norm(b); t=mul(a,conj(b)); assert t[0]%n==0 and t[1]%n==0, (a,b); return (t[0]//n,t[1]//n)
def logabs(x): return 0.5*log(norm(x))
def normalise(u):
    g=u[0]
    for x in u[1:]: g=ggcd(g,x)
    return [exdiv(x,g) for x in u]
def h_proj(u):
    u=normalise(u); return max(logabs(x) for x in u)
def h_ratio(x,a):  # Weil height of x/a in Q(i)
    g=ggcd(x,a); b=exdiv(x,g); d=exdiv(a,g); return max(logabs(b),logabs(d))
def plucker(u):
    return [x[0]*y[1]-x[1]*y[0] for x,y in combinations(u,2)]  # Im(conj(x)*y)
def h_line(u):
    p=plucker(u); G=0
    for q in p: G=gcd(G,abs(q))
    return log(max(abs(q) for q in p))-log(G), G
def minanchor(u):
    best=None
    for a in range(len(u)):
        s=sum(h_ratio(u[x],u[a]) for x in range(len(u)) if x!=a)
        if best is None or s<best[0]: best=(s,a)
    return best
def nonreal_pairs(u):
    for x,y in combinations(u,2):
        if mul(conj(x),y)[1]==0: return False
    return True
