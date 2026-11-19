# Independent sanity experiment: random parents N = product of r distinct split primes in [r^D, 2 r^D];
# enumerate ALL lattice points on x^2+y^2=N, find max number of points on an arc of length C*sqrt(R), R=sqrt(N),
# and count close pairs; compare with Poisson heuristic binom(n,2)*theta/pi, theta=C*N^(-1/4).
import random, math, sys, itertools
from sympy import isprime
def gauss_prime(p):
    # find a,b with a^2+b^2=p
    for a in range(1,int(math.isqrt(p))+1):
        b2=p-a*a; b=math.isqrt(b2)
        if b*b==b2: return (a,b)
def mul(z,w): return (z[0]*w[0]-z[1]*w[1], z[0]*w[1]+z[1]*w[0])
def run(r,D,Cs,seed):
    rnd=random.Random(seed)
    P=int(round(r**D))
    pool=[p for p in range(P,2*P+1) if p%4==1 and isprime(p)]
    if len(pool)<r: return None
    ps=rnd.sample(pool,r)
    N=1
    for p in ps: N*=p
    pts=[(1,0)]
    for p in ps:
        a,b=gauss_prime(p); pts=[mul(z,(a,b)) for z in pts]+[mul(z,(a,-b)) for z in pts]
    allp=[]
    for (x,y) in pts:
        allp += [(x,y),(-y,x),(-x,-y),(y,-x)]
    assert all(x*x+y*y==N for x,y in allp[:50])
    ang=sorted((math.atan2(y,x)%(2*math.pi),x,y) for x,y in allp)
    n=len(ang)
    res={}
    for C in Cs:
        theta=C*N**-0.25
        # exact chord test: |z_i-z_j|^2 <= C^2*R  <=> |z_i-z_j|^4 <= C^4*N  (chord vs arc: arc slightly longer;
        # use chord<=C sqrt R as a superset of arc<=C sqrt R)
        C4=round(C**4)
        best=1; pairs=0
        ext=ang+[(a+2*math.pi,x,y) for a,x,y in ang]
        j=0
        for i in range(n):
            a0,x0,y0=ext[i]
            k=i+1; cnt=1
            while k<i+n and ext[k][0]-a0 <= theta*1.001+1e-12:
                x1,y1=ext[k][1],ext[k][2]
                c2=(x1-x0)**2+(y1-y0)**2
                if c2*c2 <= C4*N:
                    cnt+=1; pairs+=1
                k+=1
            best=max(best,cnt)
        res[C]=(best,pairs,n*(n-1)/2*theta/math.pi)
    return len(pool),N,res
if __name__=="__main__":
    Cs=[1,10,100]
    for (r,D,trials) in [(10,3,20),(12,3,20),(14,3,10),(16,2.2,10),(14,4,10),(16,4,6),(18,2,4),(18,3,4)]:
        out=[]
        for t in range(trials):
            o=run(r,D,Cs,1000*r+t)
            if o is None: print(r,D,"pool too small"); break
            out.append(o)
        if not out: continue
        K=out[0][0]
        for C in Cs:
            mx=max(o[2][C][0] for o in out); tp=sum(o[2][C][1] for o in out); te=sum(o[2][C][2] for o in out)
            print(f"r={r} D={D} K={K} C={C}: max cluster over {len(out)} parents = {mx}; close pairs total={tp}, Poisson expectation={te:.3g}")
        sys.stdout.flush()
