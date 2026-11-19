# Own F_9 = F_3[i] arithmetic; elements as pairs (a,b) = a+bi; polys low-first lists
import sys, itertools
def fmul(x,y): a,b=x;c,d=y; return ((a*c-b*d)%3,(a*d+b*c)%3)
def fadd(x,y): return ((x[0]+y[0])%3,(x[1]+y[1])%3)
def fsub(x,y): return ((x[0]-y[0])%3,(x[1]-y[1])%3)
def finv(x):
    for a in range(3):
        for b in range(3):
            if fmul(x,(a,b))==(1,0): return (a,b)
def conj(x): return (x[0],(-x[1])%3)
def decode(n):
    c=[]
    while n: n,d=divmod(n,9); b,a=divmod(d,3); c.append((a,b))
    return c
def trim(p):
    p=list(p)
    while p and p[-1]==(0,0): p.pop()
    return p
def pmul(p,q):
    if not p or not q: return []
    r=[(0,0)]*(len(p)+len(q)-1)
    for i,x in enumerate(p):
        for j,y in enumerate(q): r[i+j]=fadd(r[i+j],fmul(x,y))
    return trim(r)
def psub(p,q):
    n=max(len(p),len(q)); p=p+[(0,0)]*(n-len(p)); q=q+[(0,0)]*(n-len(q))
    return trim([fsub(x,y) for x,y in zip(p,q)])
def pmod(p,q):
    p=trim(p); q=trim(q); il=finv(q[-1])
    while len(p)>=len(q):
        c=fmul(p[-1],il); s=len(p)-len(q)
        p=psub(p,[(0,0)]*s+[fmul(c,x) for x in q])
    return p
def pdivexact(p,q):
    p=trim(p); q=trim(q); il=finv(q[-1]); quo=[(0,0)]*(len(p)-len(q)+1)
    while len(p)>=len(q):
        c=fmul(p[-1],il); s=len(p)-len(q); quo[s]=c
        p=psub(p,[(0,0)]*s+[fmul(c,x) for x in q])
    assert not p; return trim(quo)
def pgcd(a,b):
    a=trim(a);b=trim(b)
    while b: a,b=b,pmod(a,b)
    il=finv(a[-1]); return [fmul(il,x) for x in a]
def deg(p): return len(trim(p))-1
def pconj(p): return [conj(c) for c in p]
def deriv(p): return trim([fmul((k%3,0),p[k]) for k in range(1,len(p))])
def show(p):
    s=[]
    for k in reversed(range(len(p))):
        a,b=p[k]
        if (a,b)==(0,0): continue
        c={(1,0):'',(2,0):'2',(0,1):'i',(0,2):'2i'}.get((a,b),f'({a}+{b}i)')
        s.append(f"{c}t^{k}" if k else (c or '1') if c else '1')
    return '+'.join(s)
# own irreducible factorisation over F_3 by trial division with monic polys (embedded in F_9)
def all_monic(d):
    for coeffs in itertools.product(range(3),repeat=d):
        yield [(c,0) for c in coeffs]+[(1,0)]
def factor_F3(N):
    N=trim(N); facs=[]; d=1
    while deg(N)>0:
        if 2*d>deg(N): facs.append((N,1)); break
        for f in all_monic(d):
            e=0
            while deg(N)>=d and not pmod(N,f):
                N=pdivexact(N,f); e+=1
            if e: facs.append((f,e))
        d+=1
    return facs
if __name__=="__main__":
  for L in open(sys.argv[1]).read().split('\n'):
    if not L.strip(): continue
    D=int(L.split()[0][2:]); M=int(L.split()[1][2:])
    marks=[decode(int(x)) for x in L.split('base9_marks=')[1].split()[0].rstrip(',').split(',')]
    assert len(marks)==M and len(set(map(tuple,marks)))==M
    norms=[pmul(z,pconj(z)) for z in marks]
    assert all(n==norms[0] for n in norms)
    N=norms[0]; assert all(c[1]==0 for c in N)
    lead=set(z[-1] for z in marks); degs=set(deg(z) for z in marks)
    gaps=[deg(psub(marks[a],marks[b])) for a,b in itertools.combinations(range(M),2)]
    sqf = deg(pgcd(N,deriv(N)))==0 if deriv(N) else False
    print(f"D={D} M={M} degs={degs} lead={lead} deg N={deg(N)} gaps={gaps} ok={2*max(gaps)<=D} squarefree={sqf}")
    print("  N =",show(N))
    if D<=16:
        fac=factor_F3(N)
        for f,e in fac:
            g=pgcd(marks[0],f) if deg(f)%2==0 else None
            if g is None or deg(g)==0 or deg(g)==deg(f):
                print(f"   prime {show(f)} (deg {deg(f)}, exp {e}) inert/unsplit-cut: gcd degs",[deg(pgcd(z,f)) for z in marks]); continue
            pat=''.join('A' if pgcd(z,f)==g else 'B' for z in marks)
            print(f"   prime {show(f)} deg {deg(f)} exp {e}: pi={show(g)}  cut {pat}")
    print("  marks:"); [print("   ",show(z)) for z in marks]
