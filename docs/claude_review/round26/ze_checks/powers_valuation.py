"""Exact check of the valuation identity behind Lemma 3.1 (forced divisibility at powers):
for z_x on x^2+y^2=N, j>=1, k in Z^n_0, |k|_1<=D, D-|k|_1 even:
  G_k(z^j) = N^{j(D-|k|_1)/2} prod (z_x^j)^{k_x^+} conj(z_x^j)^{k_x^-}  in Z[i],
  v_pi(G_k(z^j)) = j(e_p D/2 + <k,a_p>),  v_conjpi = j(e_p D/2 - <k,a_p>),
and G_k(z^j) = N^{jD/2} * Z^{jk} where Z=(z_x/z_n)."""
import itertools, random, math, cmath
def gmul(a,b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def gconj(a): return (a[0],-a[1])
def gpow(a,e):
    r=(1,0)
    for _ in range(e): r=gmul(r,a)
    return r
def vpi(x,pi,p):
    if x==(0,0): return 10**9
    v=0
    while True:
        y=gmul(x,gconj(pi))
        if y[0]%p==0 and y[1]%p==0:
            x=(y[0]//p,y[1]//p); v+=1
        else: return v
primes=[5,13,17,29]
N=math.prod(primes)
pis={}
for p in primes:
    for a in range(1,p):
        b2=p-a*a
        if b2>0 and int(math.isqrt(b2))**2==b2: pis[p]=(a,int(math.isqrt(b2))); break
pts=[(a,b) for a in range(0,int(math.isqrt(N))+1) for b in [int(math.isqrt(N-a*a))] if a*a+b*b==N and b>=0]
pts=sorted(set(pts+[(b,a) for a,b in pts]),key=lambda z: math.atan2(z[1],z[0]))
z=pts[:4]; n=len(z)
apx={p:[vpi(zz,pis[p],p) for zz in z] for p in primes}
checks=0; bad=0
random.seed(1)
for j in [1,2,3]:
    zj=[gpow(zz,j) for zz in z]
    Nj=N**j
    for trial in range(60):
        kp=[random.randint(-2,2) for _ in range(n-1)]; k=kp+[-sum(kp)]
        l1=sum(abs(t) for t in k); D=l1+2*random.randint(0,2)
        G=( N**(j*(D-l1)//2) ,0)
        for x in range(n):
            if k[x]>0: G=gmul(G,gpow(zj[x],k[x]))
            elif k[x]<0: G=gmul(G,gpow(gconj(zj[x]),-k[x]))
        for p in primes:
            e=1; ka=sum(kk*a for kk,a in zip(k,apx[p]))
            ok1 = vpi(G,pis[p],p)==j*(e*D//2+ka) if D%2==0 else True
            ok2 = vpi(G,gconj(pis[p]),p)==j*(e*D//2-ka) if D%2==0 else True
            checks+=2; bad+= (not ok1)+(not ok2)
        # archimedean identity G = N^{jD/2} Z^{jk}
        Zc=[complex(*zz)/complex(*z[-1]) for zz in z[:-1]]
        val=N**(j*D/2)*cmath.exp(sum(j*kk*cmath.log(Zc[x]) for x,kk in enumerate(kp)))
        rel=abs(complex(*G)-val)/abs(val)
        checks+=1; bad+= rel>1e-9
print("cluster",z,"valuation columns",apx)
print("checks",checks,"failures",bad)
