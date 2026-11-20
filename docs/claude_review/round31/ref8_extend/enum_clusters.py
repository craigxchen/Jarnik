# Enumerate genuine clusters: all lattice points on x^2+y^2=N for split-only N, sliding windows of M consecutive
# points (by angle), arc constant C = sqrt(R)*dtheta.  Normalise by the Gaussian gcd (subtract column minima of
# the level vectors), and record exact data for windows with C <= CMAX.
import sys, itertools, math, pickle, random
from gaussx import *
SPLIT=[p for p in range(5,400) if p%4==1 and all(p%d for d in range(2,int(p**.5)+1))]
PI={p:two_squares(p) for p in SPLIT}
def points(fac):
    ps=list(fac); es=[fac[p] for p in ps]
    pw=[[gmul(gpow(PI[p],a),gpow(gconj(PI[p]),e-a)) for a in range(e+1)] for p,e in zip(ps,es)]
    out=[]
    for a in itertools.product(*[range(e+1) for e in es]):
        z=(1,0)
        for k,ak in enumerate(a): z=gmul(z,pw[k][ak])
        for u in range(4):
            out.append((math.atan2(z[1],z[0])%(2*math.pi),z,a))
            z=(-z[1],z[0])
    out.sort()
    return ps,es,out
def windows(fac,MMAX=16,CMAX=16.0):
    ps,es,pts=points(fac); N=1
    for p,e in zip(ps,es): N*=p**e
    sR=N**0.25; n=len(pts); res=[]
    ext=pts+[(t+2*math.pi,z,a) for t,z,a in pts]
    for M in range(3,MMAX+1):
        for s in range(n):
            dth=ext[s+M-1][0]-ext[s][0]
            C=sR*dth
            if C<=CMAX: res.append((M,C,[ext[s+k][1] for k in range(M)],[ext[s+k][2] for k in range(M)]))
    return ps,es,N,res
if __name__=="__main__":
    random.seed(int(sys.argv[1]) if len(sys.argv)>1 else 1)
    NMAX=10**int(sys.argv[2]) if len(sys.argv)>2 else 10**13
    out=[]
    # norms: random products of split primes with many representations
    seen=set()
    for trial in range(int(sys.argv[3]) if len(sys.argv)>3 else 400):
        fac={}; N=1
        cand=SPLIT[:random.choice([6,8,10,14,20,30])]
        while True:
            p=random.choice(cand); 
            if N*p>NMAX: break
            N*=p; fac[p]=fac.get(p,0)+1
        key=tuple(sorted(fac.items()))
        if key in seen: continue
        seen.add(key)
        if math.prod(e+1 for e in fac.values())>6000: continue
        ps,es,N,res=windows(fac)
        out.append((ps,es,N,res))
    pickle.dump(out,open(f"clusters_{sys.argv[1] if len(sys.argv)>1 else 1}.pkl","wb"))
    print(len(out),"norms;",sum(len(r[3]) for r in out),"windows with C<=16")
