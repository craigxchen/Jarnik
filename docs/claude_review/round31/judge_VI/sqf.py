# Judge's independent check of VI.5 (squarefree binary parents), m=4 (d=3). F_2[t] as ints.
import sys, itertools
def deg(a): return a.bit_length()-1
def mul(a,b):
    r=0
    while b:
        if b&1: r^=a
        a<<=1; b>>=1
    return r
def pmod(a,m):
    dm=deg(m)
    while a and deg(a)>=dm: a^=m<<(deg(a)-dm)
    return a
def pdivmod(a,m):
    q=0; dm=deg(m)
    while a and deg(a)>=dm:
        s=deg(a)-dm; q|=1<<s; a^=m<<s
    return q,a
def gcd(a,b):
    while b: a,b=b,pmod(a,b)
    return a
def der(a):  # d/dt over F_2: keep odd powers shifted
    r=0; i=1
    while (a>>i):
        if (a>>i)&1 and i%2==1: r|=1<<(i-1)
        i+=1
    return r
def irreducibles(maxd):
    out=[]
    for d in range(1,maxd+1):
        for p in range(1<<d,1<<(d+1)):
            if all(pmod(p,q) for q in out if 2*deg(q)<=d): out.append(p)
    return out
IRR=irreducibles(10)
def factor(a):
    fs=[]
    for p in IRR:
        while a and pmod(a,p)==0 and deg(a)>=1:
            fs.append(p); a=pdivmod(a,p)[0]
        if deg(a)<1: break
    assert deg(a)<1, "factor beyond table"
    return fs
def trace(x,p):  # absolute trace of x mod p in F_2[t]/p
    e=deg(p); s=0; y=pmod(x,p)
    for _ in range(e): s^=y; y=pmod(mul(y,y),p)
    assert deg(s)<=0
    return s
d=3; m=2**(d-1); L=int(sys.argv[1]) if len(sys.argv)>1 else 16
t=2
U=[0,t,t*t,t*t^t]   # span{t,t^2}
def sq(a): return mul(a,a)
V=sorted(set(a ^ sq(sq(f)^f) for a in (0,1) for f in U))
assert len(V)==8
print("V =",[bin(v) for v in V])
Rp=sorted(set(p for v in V if v not in (0,) for p in factor(v) if deg(p)%2==0) | set())
# primes dividing a nonzero element of V (differences of elements are elements)
print("even-degree primes in R:",[bin(p) for p in Rp]," deg R =",sum(deg(p) for p in Rp))
E = (1<<L)^t if L%2==0 else (1<<L)
degR=sum(deg(p) for p in Rp)
good=[]
for S in range(1<<degR):
    T=E^sq(S)
    if all(trace(T,p)==1^((deg(p)//2)%2) for p in Rp): good.append(S)
print("number of S (deg<degR) meeting all trace conditions:",len(good),"of",1<<degR)
# F_4[t]: pairs (A,B)=A+theta B, theta^2=theta+1
def fmul(x,y):
    A,B=x; C,D=y
    AC=mul(A,C); BD=mul(B,D)
    return (AC^BD, mul(A,D)^mul(B,C)^BD)
def conj(x): A,B=x; return (A^B,B)
def fsub(x,y): return (x[0]^y[0], x[1]^y[1])
def fdeg(x): return max(deg(x[0]),deg(x[1]))
# hyperplanes of V not containing 1: kernels of functionals l with l(1)=1. basis of V: echelon
basis=[]
for v in V:
    w=v
    for b in basis:
        if w ^ b < w: w ^= b
    # simple greedy: use leading-bit elimination
basis=[]
piv={}
for v in V:
    w=v
    for b in sorted(piv, reverse=True):
        if (w>>b)&1: w^=piv[b]
    if w: piv[deg(w)]=w
basis=list(piv.values()); assert len(basis)==3
def coords(v):
    c=[0]*3; w=v
    for k,b in sorted(piv.items(), reverse=True):
        if (w>>k)&1: w^=b; c[list(piv.keys()).index(k)]=1
    assert w==0; return c
c1=coords(1)
forms=[l for l in itertools.product((0,1),repeat=3) if sum(a*b for a,b in zip(l,c1))%2==1]
assert len(forms)==4
def run(S, verbose=True):
    T=E^sq(S)
    pts=[]
    for l in forms:
        H=[v for v in V if sum(a*b for a,b in zip(l,coords(v)))%2==0]
        z=(1,0)
        for v in H:
            if v: z=fmul(z,(T^v,1))
        pts.append(z)
    Ns=[fmul(z,conj(z)) for z in pts]
    N=Ns[0]
    ok_norm=all(x==N for x in Ns) and N[1]==0
    Nn=N[0]
    sqfree = gcd(Nn,der(Nn))==1
    D=fdeg(pts[0]); assert all(fdeg(z)==D for z in pts)
    diffs=[fdeg(fsub(pts[i],pts[j])) for i in range(4) for j in range(i+1,4)]
    distinct=len(set(pts))==4
    cg=pts[0]
    # common gcd over F_4[t]: check gcd of norms of points' pairwise... simple test: gcd of the F_2-norm with conjugates done via N squarefree => no common factor with conj
    if verbose:
        print("S=",bin(S),"deg T=",deg(T)," D=",D," deg N=",deg(Nn)," norms equal & in F_2[t]:",ok_norm,
              " N squarefree:",sqfree," distinct:",distinct," pair diffs:",diffs," max<=D/2:",max(diffs)<=D/2)
    return sqfree, ok_norm, max(diffs)<=D/2
for S in good[:3]: run(S)
bad=[S for S in range(1<<degR) if S not in good]
res=[run(S,False) for S in bad]
print("control: S violating some trace condition -> squarefree in",sum(r[0] for r in res),"of",len(bad))
print("all good S squarefree & cluster:", all(all(run(S,False)) for S in good))
