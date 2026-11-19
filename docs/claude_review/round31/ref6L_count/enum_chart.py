"""Enumerate five-point clusters through the word C = A + iL of the extreme pair
(Astra's 'one chord' C, = the record's endpoint chart).  For every L in the given list and
every A <= Amax: the other three points have words X + iL with X - A | A^2 + L^2
(divisor step), pairwise compatible.  The primitive tuple is rebuilt independently by a
Gaussian lcm, and only clusters with K5 = diam / H^(1/5) <= Kmax are kept.
Usage: python3 enum_chart.py Kmax Amax L1,L2,...   (or Lmin:Lmax:step)"""
import sys, math
from math import isqrt, gcd as igcd
from fractions import Fraction
from glib import mul, conj, nrm, gcd, exdiv, lcm

def primes_upto(n):
    s=bytearray([1])*(n+1); s[0:2]=b'\0\0'
    for i in range(2,isqrt(n)+1):
        if s[i]: s[i*i::i]=bytearray(len(s[i*i::i]))
    return [i for i in range(n+1) if s[i]]

def sieve(L,Amax,primes):
    rem=[A*A+L*L for A in range(Amax+1)]
    fac=[[] for _ in range(Amax+1)]
    LL=L*L
    for p in primes:
        if p*p>rem[-1]*1.0+1 and p>Amax: break
        if p==2: roots=[L%2]
        elif L%p==0: roots=[0]
        elif p%4==3: continue
        else:
            # sqrt(-1) mod p
            c=2
            while pow(c,(p-1)//2,p)!=p-1: c+=1
            r=pow(c,(p-1)//4,p)*L%p
            roots=[r,(p-r)%p]
        for r0 in roots:
            for A in range(r0,Amax+1,p):
                if A==0: continue
                x=rem[A]; e=0
                while x%p==0: x//=p; e+=1
                if e: rem[A]=x; fac[A].append((p,e))
    for A in range(1,Amax+1):
        if rem[A]>1: fac[A].append((rem[A],1))
    return fac

def divisors_upto(f,bound):
    ds=[1]
    for p,e in f:
        nd=[]
        for d in ds:
            q=d
            for _ in range(e+1):
                if q>bound: break
                nd.append(q); q*=p
        ds=nd
    return ds

def beta_bar(X,L):
    w=(X,L); g=gcd(w,conj(w)); b=exdiv(w,g); return conj(b)

def edge_norm(X,L):
    # N(w)/N(gcd(w,conj w)), w=X+iL; gcd = g or g(1+i), g=gcd(X,L)
    g=igcd(X,L); x=X//g; l=L//g
    n=x*x+l*l
    return n//2 if (x&1 and l&1) else n

def ilcm(a,b): return a//igcd(a,b)*b

def run(L,Amax,Kmax,primes,out,seen):
    fac=sieve(L,Amax,primes); LL=L*L; ncl=0
    for A in range(1,Amax+1):
        D=A*A+LL
        T=(Kmax*math.sqrt(D)/(2*L))**(10/3)
        nA=edge_norm(A,L)
        if nA>T: continue
        Xmax=int(2*L*math.sqrt(T))+2
        cand=[]
        for d in divisors_upto(fac[A],Xmax-A):
            X=A+d
            if edge_norm(X,L)>T: continue
            Q=(A*X+LL)//d
            if edge_norm(Q,L)>T: continue
            cand.append(X)
        cand.sort()
        c=len(cand)
        base=[ilcm(nA,ilcm(edge_norm(X,L),edge_norm((A*X+LL)//(X-A),L))) for X in cand]
        adj=[dict() for _ in range(c)]
        for i in range(c):
            Xi=cand[i]; bi=base[i]
            if bi>T: continue
            for j in range(i+1,c):
                Xj=cand[j]
                q,r=divmod(Xi*Xj+LL,Xj-Xi)
                if r: continue
                v=ilcm(ilcm(bi,base[j]),edge_norm(q,L))
                if v<=T: adj[i][j]=v
        for i in range(c):
            ai=adj[i]
            for j,vij in ai.items():
                aj=adj[j]
                for k,vjk in aj.items():
                    vik=ai.get(k)
                    if vik is None: continue
                    H=ilcm(ilcm(vij,vjk),vik)
                    if H>T: continue
                    cl=[A,cand[i],cand[j],cand[k]]
                    ZG=(1,0)
                    for X in cl: ZG=lcm(ZG,beta_bar(X,L))
                    assert nrm(ZG)==H,(nrm(ZG),H,cl)
                    Delta=2*L/math.sqrt(D)
                    K5=Delta*H**0.3
                    if K5<=Kmax:
                        key=(H,frozenset(Fraction(x,L) for x in cl))
                        if key not in seen:
                            seen.add(key); ncl+=1
                            out.append((K5,H,L,tuple(cl),-math.log(Delta)/math.log(H)))
    return ncl

if __name__=='__main__':
    Kmax=float(sys.argv[1]); Amax=int(sys.argv[2]); spec=sys.argv[3]
    if ':' in spec:
        a,b,c=map(int,spec.split(':')); Ls=list(range(a,b+1,c))
    else: Ls=[int(x) for x in spec.split(',')]
    primes=primes_upto(int(math.sqrt(Amax*Amax+max(Ls)**2))+2)
    out=[];seen=set()
    for L in Ls:
        n=run(L,Amax,Kmax,primes,out,seen)
        print('L=%d new clusters %d'%(L,n),flush=True)
    out.sort()
    print('TOTAL',len(out))
    for r in out: print('K5=%.4f H=%d L=%d X=%s e=%.4f'%r)
