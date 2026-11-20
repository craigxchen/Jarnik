# For genuine windows with exponent deficiency exactly 1: unique dependence c, P=c+/|c+|, Q=c-/|c-|.
# Exact identity check: E d - W/2 = -sum_p w_p[2(mu-1/2)^2 + 2E min(1-max)] and the pair bound
# E d - W/2 >= 2E log|b| - 2 log C.  Report the three numbers.
import pickle, math, sympy
from math import log
from gaussx import *
from enum_clusters import PI
from analyse import affrank
data=[]
for f in ["scan_11_0.pkl","scan_11_1.pkl"]:
    for ps,es,N,res in pickle.load(open(f,"rb")):
        for M,C,zs,levs in res:
            cols=list(zip(*levs)); mins=[min(c) for c in cols]; maxs=[max(c) for c in cols]
            keep=[k for k in range(len(ps)) if maxs[k]>mins[k]]
            rows=[[levs[i][k]-mins[k] for k in keep] for i in range(M)]
            if M-1-affrank(rows)!=1: continue
            Mat=sympy.Matrix([[1]+r for r in rows]).T
            c=Mat.nullspace()[0]; c=[float(x) for x in c]
            sp=sum(x for x in c if x>0); P=[max(x,0)/sp for x in c]; Q=[max(-x,0)/sp for x in c]
            g=(1,0)
            for k,p in enumerate(ps): g=gmul(g,gmul(gpow(PI[p],mins[k]),gpow(gconj(PI[p]),es[k]-maxs[k])))
            zz=[gdivexact(z,g) for z in zs]
            W=sum((maxs[k]-mins[k])*log(ps[k]) for k in keep)
            Ed=0; Eb=0
            for i in range(M):
                for j in range(M):
                    if P[i]>0 and Q[j]>0:
                        d=sum(abs(rows[i][t]-rows[j][t])*log(ps[keep[t]]) for t in range(len(keep)))
                        diff=(zz[i][0]-zz[j][0],zz[i][1]-zz[j][1])
                        gij=math.gcd(*[0]) # placeholder
                        # |b|^2 = |zi-zj|^2 / |g_ij|^2, log|g_ij|^2 = W-d
                        lb=0.5*(log(gnorm(diff))-(W-d))
                        Ed+=P[i]*Q[j]*d; Eb+=P[i]*Q[j]*lb
            data.append((C,M,Ed-W/2,2*Eb-2*log(C),Eb))
print(len(data),"windows with deficiency 1")
data.sort()
for r in data[:12]: print("C=%.3f M=%d  E d-W/2=%.4f  2E log|b|-2log C=%.4f  E log|b|=%.4f"%r)
print("max of (2Elog|b|-2logC)-(Ed-W/2) (must be <=0):",max(r[3]-r[2] for r in data))
print("min E log|b| over all:",min(r[4] for r in data), " (0.5 log2 =",0.5*log(2),")")
