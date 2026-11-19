# Exact maximum cluster size on ONE parent N in F_3[t] (constant model F_9/F_3, contact floor(D/2)).
# Factor N over F_3 (sympy), split even-degree primes over F_9 by Cantor-Zassenhaus (own code),
# phases phi(pi) = pi/conj(pi) mod s^K (own code), convolve sign choices; count = #sign vectors per phase.
import sys, random, collections
from sympy import Poly, symbols
from witness import fmul,fadd,fsub,finv,conj,trim,pmul,psub,pmod,pgcd,deg,pconj,deriv
t=symbols('t')
def F3toF9(c): return [(x%3,0) for x in c]
def powmod(a,e,m):
    r=[(1,0)]; a=pmod(a,m)
    while e:
        if e&1: r=pmod(pmul(r,a),m)
        a=pmod(pmul(a,a),m); e>>=1
    return r
def split(p):  # p in F_9[t] monic product of two conjugate irreducibles of degree e
    e=deg(p)//2
    while True:
        a=[(random.randrange(3),random.randrange(3)) for _ in range(deg(p))]
        a=trim(a)
        if deg(a)<1: continue
        g=pgcd(a,p)
        if 0<deg(g)<deg(p): f=g
        else:
            h=psub(powmod(a,(9**e-1)//2,p),[(1,0)])
            if not h: continue
            f=pgcd(h,p)
        if deg(f)==e: return f
def factor_parent(Nlist):  # Nlist low-first ints
    P=Poly(list(reversed(Nlist)),t,modulus=3)
    lc,fl=P.factor_list()
    primes=[]
    for f,ex in fl:
        if ex!=1: return None
        c=[int(x)%3 for x in reversed(f.all_coeffs())]
        if (len(c)-1)%2: return None   # inert prime: odd valuation impossible for a norm
        primes.append(split(F3toF9(c)))
    return primes
def phase(p,K):
    d=deg(p); a=[p[d-i] if i<=d else (0,0) for i in range(K)]; b=[conj(x) for x in a]
    q=[]
    for i in range(K):
        v=a[i]
        for j in range(1,i+1): v=fsub(v,fmul(b[j],q[i-j]))
        q.append(v)
    return tuple(q)
def smul(x,y,K):
    r=[(0,0)]*K
    for i in range(K):
        if x[i]==(0,0): continue
        for j in range(K-i): r[i+j]=fadd(r[i+j],fmul(x[i],y[j]))
    return tuple(r)
def sinv(x,K):
    q=[]
    for i in range(K):
        v=(1,0) if i==0 else (0,0)
        for j in range(1,i+1): v=fsub(v,fmul(x[j],q[i-j]))
        q.append(v)
    return tuple(q)
def parent_max(Nlist,D):
    K=(D+1)//2
    primes=factor_parent(Nlist)
    if primes is None: return None
    cnt=collections.Counter({tuple([(1,0)]+[(0,0)]*(K-1)):1})
    for p in primes:
        f=phase(p,K); fi=sinv(f,K); new=collections.Counter()
        for x,c in cnt.items(): new[smul(x,f,K)]+=c; new[smul(x,fi,K)]+=c
        cnt=new
    m=max(cnt.values()); return m, sum(1 for v in cnt.values() if v==m), len(primes)
if __name__=="__main__":
    N=[int(x) for x in sys.argv[1].split(',')]; D=(len(N)-1)//2
    print(parent_max(N,D))

import subprocess
def parent_max_fast(Nlist,D):
    K=(D+1)//2
    primes=factor_parent(Nlist)
    if primes is None: return None
    lines=[f"{K} {len(primes)}"]
    for p in primes:
        f=phase(p,K); fi=sinv(f,K)
        lines.append(' '.join(str(a+3*b) for a,b in f)+' '+' '.join(str(a+3*b) for a,b in fi))
    out=subprocess.run(['./conv'],input='\n'.join(lines),capture_output=True,text=True).stdout.split()
    return int(out[0]),int(out[1]),len(primes)
