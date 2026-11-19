# Independent check of Part VI.5 (squarefree binary parents). Own F2[t]/F4[t] arithmetic.
# F2[t]: python int, bit i = coeff of t^i.  F4[t]: pair (A,B) meaning A + B*theta, theta^2=theta+1.
import sys, random, itertools
def deg(a): return a.bit_length()-1
def mul(a,b):
    if a.bit_length()<b.bit_length(): a,b=b,a
    r=0
    while b:
        low=b&-b; r^=a<<(low.bit_length()-1); b^=low
    return r
def mod(a,p):
    dp=deg(p)
    while a and deg(a)>=dp: a^=p<<(deg(a)-dp)
    return a
def gcd(a,b):
    while b: a,b=b,mod(a,b)
    return a
def der(a):
    r=0;i=1
    x=a>>1; # coefficient of t^i (i odd) -> t^(i-1)
    mask=0
    n=a.bit_length()
    for j in range(0,n,2): mask|=1<<j
    return x&mask
def sq(a): return mul(a,a)
def powmod(a,e,p):
    r=1;a=mod(a,p)
    while e:
        if e&1: r=mod(mul(r,a),p)
        a=mod(mul(a,a),p); e>>=1
    return r
def irreducible(p):
    n=deg(p)
    if n<=0: return False
    # Rabin: t^(2^n)=t mod p and gcd(t^(2^(n/q))-t,p)=1 for prime q|n
    x=2
    def frob(k):
        y=2
        for _ in range(k): y=mod(mul(y,y),p)
        return y
    if frob(n)!=mod(2,p): return False
    for q in [q for q in range(2,n+1) if n%q==0 and all(q%s for s in range(2,q))]:
        if gcd(p, frob(n//q)^2)!=1: return False
    return True
def trace(a,p):
    a=mod(a,p); s=0
    for _ in range(deg(p)): s^=a; a=mod(mul(a,a),p)
    assert s in (0,1); return s
# F4[t]
def fmul(x,y):
    A1,B1=x;A2,B2=y
    bb=mul(B1,B2)
    return (mul(A1,A2)^bb, mul(A1,B2)^mul(A2,B1)^bb)
def fconj(x): A,B=x; return (A^B,B)
def fdeg(x): return max(deg(x[0]),deg(x[1]))
def fprod(lst):
    # balanced product
    lst=list(lst)
    if not lst: return (1,0)
    while len(lst)>1:
        lst=[fmul(lst[i],lst[i+1]) if i+1<len(lst) else lst[i] for i in range(0,len(lst),2)]
    return lst[0]
def prod2(lst):
    lst=list(lst)
    while len(lst)>1:
        lst=[mul(lst[i],lst[i+1]) if i+1<len(lst) else lst[i] for i in range(0,len(lst),2)]
    return lst[0]

def build_V(d):
    # basis: 1, (f^2+f)^2 for f=t^1..t^(d-1); V = all F2-combinations (map additive)
    basis=[1]+[sq(sq(1<<j)^(1<<j)) for j in range(1,d)]
    V=[]
    for c in range(1<<d):
        v=0
        for i in range(d):
            if c>>i&1: v^=basis[i]
        V.append(v)
    assert len(set(V))==1<<d
    # check each v = a+(f^2+f)^2 is a square: derivative 0
    assert all(der(v)==0 for v in V)
    return basis,V

def even_primes_R(d,V):
    maxdeg=2*(d-1)
    ps=[]
    for p in range(4,1<<(maxdeg+1)):
        if deg(p)%2==0 and irreducible(p) and any(v and mod(v,p)==0 for v in V):
            ps.append(p)
    return ps

def check(d,L=None,seed=0,trace_cond=True,verbose=True):
    m=1<<(d-1)
    basis,V=build_V(d)
    ps=even_primes_R(d,V)
    R=prod2(ps) if ps else 1
    # trace zero claim
    for p in ps:
        for v in V: assert trace(v,p)==0
    # check every even-degree prime of degree<=2(d-1) dividing a nonzero v (all such are in ps by construction)
    # H bound
    H=prod2([sq(f)^f^1 for f in range(0,1<<d,2)])  # f in span(t..t^(d-1)) = even ints < 2^d
    assert deg(H)==2+2*m*(d-2)
    assert mod(H,R)==0
    if L is None: L=4*m*(d-2)
    # choose S by CRT, independently: solve per prime by brute search over residues, then combine
    rnd=random.Random(seed)
    E=(1<<L)^2 if L%2==0 else (1<<L)
    S=0;Mod=1
    for p in ps:
        e=deg(p); trth=(e//2)%2
        goal=1^trth
        # want Tr(E+S^2)=goal mod p; Tr(S^2)=Tr(S)
        need=goal^trace(E,p)
        # find residue s0 with Tr(s0)=need (random)
        while True:
            s0=rnd.randrange(1<<e)
            if trace(s0,p)==need: break
        # CRT: S'=S+Mod*k with S' = s0 mod p
        inv=powmod(Mod,(1<<e)-2,p)
        k=mod(mul(s0^mod(S,p),inv),p)
        S^=mul(Mod,k); Mod=mul(Mod,p)
        assert mod(S,p)==s0
    if not trace_cond:
        S=rnd.randrange(1<<max(1,deg(R)))
    T=E^sq(S)
    assert deg(T)==L
    assert L>max(deg(v) for v in V)
    theta=(0,1)
    fac={v:(T^v,1) for v in V}   # T+theta+v
    N=prod2([sq(T^v)^(T^v)^1 for v in V if v not in (0,1) and not (v&1)])  # pairs {v,v+1}; norm of (A+theta)=A^2+A+1
    NF=fprod([fac[v] for v in V if v not in (0,1)]); assert NF==(N,0)
    sqfree = gcd(N,der(N))==1
    res=dict(d=d,m=m,L=L,degR=deg(R),nprimes=len(ps),bound=2*m*(d-2),degN=deg(N),sqfree=sqfree)
    if not trace_cond: return res
    assert sqfree
    # marks: forms l with l(1)=1: l given by bitmask over basis coords, bit0 set
    coords={v:c for c,v in enumerate(V)}
    forms=[l for l in range(1<<d) if l&1]
    marks=[]
    for l in forms:
        Hs=[v for v in V if bin(coords[v]&l).count('1')%2==0]
        z=fprod([fac[v] for v in Hs if v!=0])
        marks.append((l,Hs,z))
    D=(m-1)*L
    for l,Hs,z in marks:
        assert fdeg(z)==D
        A,B=z
        n=sq(A)^mul(A,B)^sq(B)
        assert n==N
    assert len(set(z for _,_,z in marks))==m
    maxgap=0;npairs=0
    for (l1,H1,z1),(l2,H2,z2) in itertools.combinations(marks,2):
        g=fdeg((z1[0]^z2[0],z1[1]^z2[1]))
        J=set(H1)&set(H2)
        PJ1=prod2([1^j for j in J])
        assert g==(m//2-1)*L+deg(PJ1)
        assert 2*g<=D
        maxgap=max(maxgap,g);npairs+=1
    # Walsh-Hadamard: E0 = span of basis[1:], blocks v in E0\{0}; mark l chooses v (if v in H) else v+1
    E0=[v for v in V if coords[v]&1==0]
    rows=[]
    for l,Hs,z in marks:
        Hset=set(Hs)
        row=[]
        for v in E0:
            if v==0: continue
            assert (v in Hset)!=((v^1) in Hset)
            row.append(1 if v in Hset else -1)
        rows.append([1]+row)
    for i in range(m):
        for j in range(m):
            s=sum(a*b for a,b in zip(rows[i],rows[j]))
            assert s==(m if i==j else 0)
    # disjoint prime supports of blocks N_v=f_v*fbar_v
    Nv={v:mul(T^v,T^v)^(T^v)^1 for v in E0 if v}
    for a,b in itertools.combinations(list(Nv),2):
        assert gcd(Nv[a],Nv[b])==1
    for v in Nv:  # f_v coprime to conjugate: N_v squarefree
        assert gcd(Nv[v],der(Nv[v]))==1
    res.update(D=D,maxgap=maxgap,pairs=npairs,maxPJ1=maxgap-(m//2-1)*L)
    return res

if __name__=="__main__":
    for d in map(int,sys.argv[1:]):
        print(check(d),flush=True)
