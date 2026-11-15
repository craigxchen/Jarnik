# Independent exact position-average of E^W: A = sum_n sum_{z!=w, |z|^2=|w|^2=n} (Theta-angle(z,w))_+/(2 pi)
# Enumerate first-octant-free: all points in annulus via vectorised numpy; group by norm with argsort.
import sys, math, numpy as np
def A_of(N, C, W):
    R = math.sqrt(N); Th = C*math.sqrt(R)/R
    r = math.isqrt(N+W)+1
    xs=[];ys=[]
    for x in range(-r, r+1):
        lo=N-W-x*x; hi=N+W-x*x
        if hi<0: continue
        ylo=0 if lo<=0 else math.isqrt(lo-1)+1; yhi=math.isqrt(hi)
        y=np.arange(ylo,yhi+1)
        xs.append(np.full(len(y),x)); ys.append(y)
        y2=y[y>0]; xs.append(np.full(len(y2),x)); ys.append(-y2)
    x=np.concatenate(xs).astype(np.int64); y=np.concatenate(ys).astype(np.int64)
    n=x*x+y*y; ang=np.arctan2(y,x)
    o=np.lexsort((ang,n)); n=n[o]; ang=ang[o]
    A=0.0
    # pairs within same norm: loop over offset d while any same-norm pair within Theta
    L=len(n)
    for d in range(1,200):
        same=(n[d:]==n[:-d])
        if not same.any(): break
        diff=(ang[d:]-ang[:-d])[same]   # in (0,2pi), sorted so positive
        sep=np.minimum(diff,2*math.pi-diff)
        A+=2*np.clip(Th-sep,0,None).sum()/(2*math.pi)
    return A, len(n)
if __name__=="__main__":
    for N in map(int, sys.argv[1:]):
        R=math.sqrt(N)
        for C in (4,8):
            for wf in (1.0,0.1):
                W=int(wf*R); A,npts=A_of(N,C,W)
                base=(W/R)*C*C
                print(f"N={N} logR={math.log(R):.3f} C={C} W/R={W/R:.3f} pts={npts} A={A:.3f} A/((W/R)C^2 logR)={A/(base*math.log(R)):.4f} A/((W/R)C^2)={A/base:.4f} heur(2/pi^2)(logR/2+logC)={2/math.pi**2*(math.log(R)/2+math.log(C)):.4f}",flush=True)
