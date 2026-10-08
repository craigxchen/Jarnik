import sys, math
sys.path.insert(0,'.')
from math import isqrt
import ref_check as rc
def factor(n):
    f={};d=2
    while d*d<=n:
        while n%d==0: f[d]=f.get(d,0)+1; n//=d
        d+=1
    if n>1: f[n]=f.get(n,0)+1
    return f
def best_cluster(n,k):
    P=[(x,isqrt(n-x*x)) for x in range(1,isqrt(n)+1) if isqrt(n-x*x)**2==n-x*x]
    A=sorted((math.atan2(y,x),(x,y)) for x,y in P)
    c=len(A); best=None
    for s in range(c):
        pts=[]
        for j in range(k):
            e=s+j; th,p=A[e%c]; q=p
            for _ in range(e//c): q=(-q[1],q[0])
            pts.append(q)
        C,N,kk,gn,Ps=rc.cluster_C(pts,reduce_gcd=False)
        Cr=rc.cluster_C(pts)[0]
        if best is None or C<best[0]: best=(C,Cr,pts,gn)
    return best
for n,k in [(int(a),int(b)) for a,b in (s.split(':') for s in sys.argv[1:])]:
    C,Cr,pts,gn=best_cluster(n,k)
    print("n=%d %s k=%d rawC=%.6f reducedC=%.6f gcdnorm=%d pts=%s"%(n,factor(n),k,C,Cr,gn,pts))
