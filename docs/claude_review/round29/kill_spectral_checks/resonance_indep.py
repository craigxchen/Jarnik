# Exact E^W(I) at the predicted resonant arc and sup over all positions (sweep), W=5R and W=R.
import sys, math, numpy as np
def pts(N,W):
    r=math.isqrt(N+W)+1; xs=[];ys=[]
    for x in range(-r,r+1):
        lo=N-W-x*x; hi=N+W-x*x
        if hi<0: continue
        ylo=0 if lo<=0 else math.isqrt(lo-1)+1; yhi=math.isqrt(hi)
        y=np.arange(ylo,yhi+1); xs.append(np.full(len(y),x)); ys.append(y)
        y2=y[y>0]; xs.append(np.full(len(y2),x)); ys.append(-y2)
    x=np.concatenate(xs).astype(np.int64); y=np.concatenate(ys).astype(np.int64)
    return x,y
def run(N,C,wf):
    R=math.sqrt(N); W=int(wf*R); Th=C*math.sqrt(R)/R
    x,y=pts(N,W); n=x*x+y*y; ang=np.mod(np.arctan2(y,x),2*math.pi)
    o=np.lexsort((ang,n)); n=n[o]; ang=ang[o]
    starts=[];ends=[]
    for d in range(1,400):
        same=(n[d:]==n[:-d])
        if not same.any(): break
        a1=ang[:-d][same]; a2=ang[d:][same]   # a1<a2
        for lo,hi in ((a1,a2),(a2,a1+2*math.pi)):
            ok=(hi-lo)<Th
            starts.append(hi[ok]-Th); ends.append(lo[ok])
    s=np.mod(np.concatenate(starts),2*math.pi); e=s+(np.concatenate(ends)-(np.concatenate(starts)))
    # interval [s, e] (length Th - sep), weight 2 (ordered pair); duplicate shifted by 2pi for wrap
    ev_t=np.concatenate([s,e,s+2*math.pi,e+2*math.pi]); ev_w=np.concatenate([2*np.ones(len(s)),-2*np.ones(len(s))]*2)
    o=np.lexsort((-ev_w,ev_t)); cum=np.cumsum(ev_w[o]); t=ev_t[o]
    keep=(t>=2*math.pi)&(t<4*math.pi)
    tt=np.mod(t[keep],2*math.pi); cc=cum[keep]
    # nontrivial: distance >= 2Th from multiples of pi/4 for the whole arc [theta,theta+Th]
    def nontriv(th):
        dd=np.minimum(np.mod(th,math.pi/4),math.pi/4-np.mod(th,math.pi/4))
        dd2=np.minimum(np.mod(th+Th,math.pi/4),math.pi/4-np.mod(th+Th,math.pi/4))
        return (dd>=2*Th)&(dd2>=2*Th)&((np.mod(th,math.pi/4)+Th)<math.pi/4)
    nt=nontriv(tt)
    sup_all=int(cc.max()); sup_nt=int(cc[nt].max()); th_star=tt[nt][np.argmax(cc[nt])]
    # predicted resonant arc centred at arg(i(1+2i)) = pi - atan(1/2)
    c0=math.pi-math.atan(0.5); th0=c0-Th/2
    sep_in=lambda a,b: (np.mod(a-th0,2*math.pi)<=Th)&(np.mod(b-th0,2*math.pi)<=Th)
    # direct E at th0
    inarc=np.mod(ang-th0,2*math.pi)<=Th; nn=n[inarc]
    u,cnt=np.unique(nn,return_counts=True); E0=int((cnt*(cnt-1)).sum())
    pred=2*math.floor(min(C/(2*math.sqrt(5)),1)*math.sqrt(R))
    print(f"N={N} C={C} W/R={wf} sqrtR={math.sqrt(R):.1f} E(pred arc)={E0} lowerbound_pred={pred} E/sqrtR={E0/math.sqrt(R):.3f} sup_nontriv={sup_nt} at {math.degrees(th_star):.2f}deg sup_all={sup_all}",flush=True)
for N in map(int,sys.argv[1:]):
    for wf in (5.0,1.0): run(N,4,wf)
