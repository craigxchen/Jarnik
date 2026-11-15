# Kill of Prop 4.4: inside one class mod q, residues mod l^(t+b) are confined to l^b classes,
# so the depth-one pigeonhole of two_thirds/proof.md Sec.6 survives at every l | q.
import math, sympy, kc
from math import gcd
from fractions import Fraction
def E(n,Q):
    b,r=divmod(n,Q); return Q*b*(b-1)//2+r*b
# (A) exact lifting count
bad=0; tot=0
for l in [3,5,7,11,13]:
    for N in [1,2,5,13,17,29,65,85,221,1105]:
        if N%l==0: continue
        for t in [1,2]:
            for b in [1,2]:
                if l**(t+b)>3000: continue
                m=l**t; MM=l**(t+b)
                for x1 in range(m):
                    for y1 in range(m):
                        if (x1*x1+y1*y1-N)%m: continue
                        c=sum(1 for x in range(x1,MM,m) for y in range(y1,MM,m) if (x*x+y*y-N)%MM==0)
                        tot+=1; bad+= (c!=l**b)
print("(A) lifting count = l^b:",tot,"base classes, failures",bad)
# Gaussian helpers
def gdiv_exact(a,b):
    n=b[0]**2+b[1]**2; re=a[0]*b[0]+a[1]*b[1]; im=a[1]*b[0]-a[0]*b[1]
    assert re%n==0 and im%n==0; return (re//n,im//n)
def gmod(a,b):
    n=b[0]**2+b[1]**2; re=a[0]*b[0]+a[1]*b[1]; im=a[1]*b[0]-a[0]*b[1]
    qr=(2*re+n)//(2*n); qi=(2*im+n)//(2*n); p=kc.gmul(b,(qr,qi)); return (a[0]-p[0],a[1]-p[1])
def ggcd(a,b):
    while b!=(0,0): a,b=b,gmod(a,b)
    return a
def nrm(z): return z[0]**2+z[1]**2
def vpi(z,pi):
    v=0
    while True:
        n=nrm(pi); re=z[0]*pi[0]+z[1]*pi[1]; im=z[1]*pi[0]-z[0]*pi[1]
        if re%n or im%n: return v
        z=(re//n,im//n); v+=1
# (B) smooth-modulus classes
import random
split=[p for p in sympy.primerange(5,400) if p%4==1]
random.seed(7)
Ns=[math.prod(split[:k]) for k in range(3,9)]+[5**3*13**2*17**2,5**2*13**2*17**2*29]
for _ in range(8):
    Ns.append(math.prod(p**random.choice([1,1,2]) for p in random.sample(split[:25],random.randint(4,7))))
cnt=0; failA=0; failB=0; failC=0; failD=0; maxM=0
for N in Ns:
    P=kc.points(N)
    if len(P)>2500: continue
    for lam in [1.0,0.5,0.25]:
        Q=lam*N**0.25; q=1
        for l in sympy.primerange(3,10**5):
            if N%l: q*=l
            if q>=Q: break
        fac=sympy.factorint(q)
        for key,cl in kc.classes(P,q).items():
            M=len(cl)
            if M<2: continue
            cnt+=1; maxM=max(maxM,M); K=M*(M-1)//2
            for l,t in fac.items():
                res=set((x%l**(t+1),y%l**(t+1)) for x,y in cl)
                if len(res)>l: failA+=1
                s=sum(sympy.multiplicity(l,gcd(zi[0]-zj[0],zi[1]-zj[1]))-t for i,zi in enumerate(cl) for zj in cl[i+1:])
                if s<E(M,l): failB+=1
            # normalised cut inequality, Theorem 4.3(i)
            g=cl[0]
            for z in cl[1:]: g=ggcd(g,z)
            zs=[gdiv_exact(z,g) for z in cl]; Np=N//nrm(g); Rp=math.sqrt(Np)
            lhs=0.0
            for i in range(M):
                for j in range(i+1,M):
                    h=ggcd(zs[i],zs[j]); c=gdiv_exact((zs[i][0]-zs[j][0],zs[i][1]-zs[j][1]),h)
                    assert c[0]%q==0 and c[1]%q==0
                    lhs+=0.5*math.log(nrm(c))-math.log(q)
            Qt=0.0
            for p,e in sympy.factorint(Np).items():
                pi=kc.gprime(p); a=[vpi(z,pi) for z in zs]
                for tau in range(1,e+1):
                    S=sum(1 for x in a if x<tau); Qt+=math.log(p)*(S-M/2)**2
            lamp=q/Np**0.25
            if lhs+Qt/2> M/4*math.log(Rp)+K*math.log(2/lamp)+1e-9: failC+=1
            if lhs< K*math.log(2)/2 -1e-9: failD+=1
print("(B) classes with >=2 points:",cnt,"max size",maxM,"| confinement failures",failA,"| pigeonhole E(M,l) failures",failB,"| Thm4.3(i) failures",failC,"| ramified lower bound failures",failD)
