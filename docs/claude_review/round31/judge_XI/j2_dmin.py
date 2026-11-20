import math,sys
def isprime(n): return n>1 and all(n%d for d in range(2,math.isqrt(n)+1))
def two_sq(p):
    for a in range(1,math.isqrt(p)+1):
        b=math.isqrt(p-a*a)
        if a*a+b*b==p: return (a,b)
def gm(x,y): return (x[0]*y[0]-x[1]*y[1], x[0]*y[1]+x[1]*y[0])
for P,res in [(200,1),(200,5)]:
    primes=[p for p in range(P,2*P+1) if isprime(p) and (p%4==1 if res==1 else p%8==5)]
    q=math.prod(primes)
    pts=[(1,0)]
    for p in primes:
        a,b=two_sq(p); pts=[gm(z,f) for z in pts for f in ((a,b),(a,-b))]
    def rot(z):
        x,y=z
        while not (x>0 and y>=0): x,y=-y,x
        return (x,y)
    pts=[rot(z) for z in pts]
    K=4*q.bit_length()+20
    pts.sort(key=lambda z:(z[1]<<K)//z[0])
    for a,b in zip(pts,pts[1:]): assert a[0]*b[1]-a[1]*b[0]>0
    pts2=pts+[(-pts[0][1],pts[0][0])]
    L=min((a[0]-b[0])**2+(a[1]-b[1])**2 for a,b in zip(pts2,pts2[1:]))
    fl=math.isqrt(math.isqrt(L*L//q))  # floor((L^2/q)^{1/4}) approx exact floor of dmin/q^{1/4}
    # exact floor: largest k with k^4 q <= L^2
    k=fl
    while (k+1)**4*q<=L*L: k+=1
    while k**4*q>L*L: k-=1
    print(f"P={P} family={'1mod4' if res==1 else '5mod8'} s={len(primes)} dmin^2={L} floor(dmin/q^1/4)={k} M_max(q,8)=1 certified:{L*L>8**4*q}")
