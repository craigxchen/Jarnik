# Judge's own check of the cocycle referee's q_200 triple (R5), class counts (XI.5), rigidity, pinning.
import math, random, itertools
from fractions import Fraction
def isprime(n):
    if n<2: return False
    for d in range(2,int(n**0.5)+1):
        if n%d==0: return False
    return True
def two_sq(p):
    for a in range(1,int(p**0.5)+1):
        b2=p-a*a; b=math.isqrt(b2)
        if b*b==b2: return (a,b)
def gmul(x,y): return (x[0]*y[0]-x[1]*y[1], x[0]*y[1]+x[1]*y[0])
def conj(x): return (x[0],-x[1])
P=200
primes=[p for p in range(P,2*P+1) if isprime(p) and p%4==1]
q=math.prod(primes); s=len(primes)
print("s=",s,"q=",q)
pis=[two_sq(p) for p in primes]
# all norm-q points (one representative per unit orbit, then 4 units)
pts=[(1,0)]
for pi in pis:
    pts=[gmul(z,f) for z in pts for f in (pi,conj(pi))]
units=[(1,0),(0,1),(-1,0),(0,-1)]
allpts=set(gmul(u,z) for z in pts for u in units)
assert len(allpts)==4*2**s and all(x*x+y*y==q for x,y in allpts)
G=4391525450
z=[(0,0),(3439204080,29559507552),(11203893610,21049952704)]
found=[]
for A in allpts:
    B=(A[0]+G*z[1][0],A[1]+G*z[1][1]); Cc=(A[0]+G*z[2][0],A[1]+G*z[2][1])
    if B in allpts and Cc in allpts: found.append((A,B,Cc))
print("canonical triples A_i = A + G z_i on the circle:",len(found))
def row(A):
    r=[]
    for p,pi in zip(primes,pis):
        t=gmul(A,conj(pi))
        if t[0]%p==0 and t[1]%p==0: r.append(1)
        else:
            t2=gmul(A,pi); assert t2[0]%p==0 and t2[1]%p==0; r.append(0)
    return r
def unit_of(A,r):
    prod=(1,0)
    for b,pi in zip(r,pis): prod=gmul(prod,pi if b else conj(pi))
    # A = eta*prod
    for u in units:
        if gmul(u,prod)==A: return u
qq4=q  # use integer test 16 t^4 q <= C^4 Q^2
for (A1,A2,A3) in found[:1]:
    As=[A1,A2,A3]; rows=[row(a) for a in As]; etas=[unit_of(a,r) for a,r in zip(As,rows)]
    print("units:",etas)
    Cmin=0
    for i,j in itertools.combinations(range(3),2):
        Q=math.prod(p for p,a,b in zip(primes,rows[i],rows[j]) if a!=b)
        H=(1,0)
        for p,pi,a,b in zip(primes,pis,rows[i],rows[j]):
            if a==b: H=gmul(H,pi if a else conj(pi))
        d=(As[j][0]-As[i][0],As[j][1]-As[i][1])
        den=gmul(gmul((0,2),etas[i]),H)   # 2 i eta H
        nH=den[0]**2+den[1]**2
        num=gmul(d,conj(den))
        assert num[1]==0 and num[0]%nH==0
        T=num[0]//nH
        dz=(z[j][0]-z[i][0],z[j][1]-z[i][1])
        numz=gmul(dz,conj(den)); assert numz[1]==0 and numz[0]%nH==0
        t=numz[0]//nH
        assert t*G==T
        ratio=2*abs(t)*q**0.25/math.sqrt(Q)
        Cmin=max(Cmin,ratio)
        # pinning: canonical squared chord residue kept at p in D_ij?
        kept=[p for p,a,b in zip(primes,rows[i],rows[j]) if a!=b and (G*G-1)%p==0]
        print(f"pair {i}{j}: #D={sum(a!=b for a,b in zip(rows[i],rows[j]))} T={T} t={t} ratio={ratio:.4f} pass C=4:{16*t**4*q<=4**4*Q*Q} C=3:{16*t**4*q<=3**4*Q*Q} kept-canonical-residue primes:{kept}")
    print("minimal C = %.4f"%Cmin)
    print("gcd(G^2-1, q) =",math.gcd(G*G-1,q))
# class counts mod l^(a+1) for norm-q points sharing a class mod l^a (XI.5), exact brute force
for l in [3,5,7,11,13]:
    for a in [0,1,2]:
        mod=l**(a+1); low=l**a
        for qq in [x for x in [q, 5*13*17, 29*37, 41] if x%l]:
            sols=[(x,y) for x in range(mod) for y in range(mod) if (x*x+y*y-qq)%mod==0]
            if a==0:
                assert len(sols)==l-(1 if l%4==1 else -1), (l,a,len(sols))
            else:
                from collections import Counter
                c=Counter((x%low,y%low) for x,y in sols)
                assert set(c.values())=={l}, (l,a,set(c.values()))
print("class counts: l-(-1/l) at a=0, exactly l lifts at a=1,2: OK for l=3..13")
