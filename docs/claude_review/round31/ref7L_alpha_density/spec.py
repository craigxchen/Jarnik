# Integer specialisation with exact Gaussian-integer arithmetic (pairs of Python ints).
from fractions import Fraction as F
import math, sys
EDGES = [(0,1),(0,2),(1,2),(0,3)]
def mul(a,b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def cj(a): return (a[0],-a[1])
def nrm(a): return a[0]*a[0]+a[1]*a[1]
def sub(a,b): return (a[0]-b[0],a[1]-b[1])
def gdivmod(a,b):
    n=nrm(b); num=mul(a,cj(b))
    q=(round(F(num[0],n)), round(F(num[1],n)))
    return q, sub(a,mul(q,b))
def ggcd(a,b):
    while b!=(0,0):
        _,r=gdivmod(a,b); a,b=b,r
    return a
def gexact(a,b):
    q,r=gdivmod(a,b); assert r==(0,0); return q
def matpow(T,n):
    A=[[T*T+1,T],[T,1]]; M=[[1,0],[0,1]]
    for _ in range(n):
        M=[[M[0][0]*A[0][0]+M[0][1]*A[1][0], M[0][0]*A[0][1]+M[0][1]*A[1][1]],
           [M[1][0]*A[0][0]+M[1][1]*A[1][0], M[1][0]*A[0][1]+M[1][1]*A[1][1]]]
    return M
def tuple_at(n,seeds,T):
    M=matpow(T,n); G=[]
    for s in seeds:
        a,b=s.numerator,s.denominator   # cleared: (b, a) instead of (1, a/b)
        X=M[0][0]*b+M[0][1]*a; Y=M[1][0]*b+M[1][1]*a
        G.append((X,Y))
    beta=(T,2); z0=beta
    for g in G: z0=mul(z0,cj(g))
    zs=[z0]
    for (j,k) in EDGES:
        z=mul(cj(beta),mul(G[j],G[k]))
        for l in range(4):
            if l not in (j,k): z=mul(z,cj(G[l]))
        zs.append(z)
    return zs
def analyse(n,seeds,T):
    zs=tuple_at(n,seeds,T)
    N=nrm(zs[0]); assert all(nrm(z)==N for z in zs)
    if len(set(zs))<5: return None
    g=zs[0]
    for z in zs[1:]: g=ggcd(g,z)
    red=[gexact(z,g) for z in zs]
    H=nrm(red[0])                     # H_{-K} = R''^2 of reduced tuple
    assert H*nrm(g)==N
    d2=max(F(nrm(sub(red[j],red[0])),H) for j in range(1,5))   # dist^2 ~ max|zeta_j/zeta_1-1|^2
    gam=F(8*n+1,2*n)
    logval = float(gam/2)*(math.log(d2.numerator)-math.log(d2.denominator)) + math.log(H)   # log(dist^gamma * H)
    expo = -math.log(H)/(0.5*(math.log(d2.numerator)-math.log(d2.denominator)))           # gamma with dist^gamma H = 1
    return nrm(g), H, logval, expo
if __name__=='__main__':
    seedsets=[(F(1,2),F(-2),F(7,3),F(4)),(F(-1,3),F(5,2),F(6),F(-7)),(F(0),F(1),F(3,5),F(-11,4))]
    for n in (1,2,3):
        for sd in seedsets:
            gs=set(); rows=[]
            for T in list(range(2,200))+[10**3,10**4,10**6,10**9,10**15]:
                r=analyse(n,sd,T)
                if r is None: print('   skip non-distinct T=',T); continue
                gn,H,lv,ex=r; gs.add(gn)
                if T in (10,100,10**3,10**6,10**9,10**15): rows.append((T,round(lv,4),round(ex,5)))
            print('n',n,'seeds',[str(x) for x in sd],'N(gcd) values',sorted(gs)[:12], 'gamma0=',float(F(8*n+1,2*n)))
            for r in rows: print('   T=%g  log(dist^g0*H)=%s  gamma_eff=%s'%r)
