# zero-sum sets of Q_cube(h(t)) at D=16, c=0 (L=7); numerical enumeration of all 2^16 subsets, then exact check
import numpy as np, itertools, mpmath as mp, sys
mp.mp.dps=40
y=[mp.mpf(-9)/4,mp.mpf(-2)/3,mp.mpf(1)/4,mp.mpf(8)/3]; eps=[1,-1,-1,1]
r=[eps[k]/y[k]+1j*y[k] for k in range(4)]; r=[mp.mpc(x) for x in r]
base=r+[-x for x in r]
def roots_of_comp(a,b):
    # h(t)=t^2+a t + b ; solve h(t)=rho
    out=[]
    for rho in base:
        disc=mp.sqrt(a*a-4*(b-rho))
        out+= [(-a+disc)/2,(-a-disc)/2]
    return out
L=7
for (a,b) in [(0,0),(1,0),(mp.mpf(3)/7,mp.mpf(-2)/5)]:
    rts=roots_of_comp(a,b)
    V=np.array([[float(mp.im(x**m)) for x in rts] for m in range(1,L+1)])  # L x 16
    scale=np.abs(V).max(axis=1,keepdims=True); V=V/scale
    D=16
    subs=np.array(list(itertools.product([0,1],repeat=D)),dtype=float)
    S=subs@V.T
    ok=np.nonzero(np.abs(S).max(axis=1)<1e-8)[0]
    # exact-ish recheck at 40 digits
    good=[]
    for i in ok:
        F=[k for k in range(D) if subs[i,k]]
        if all(abs(sum(mp.im(rts[k]**m) for k in F))<mp.mpf(10)**-25 for m in range(1,L+1)): good.append(F)
    rk=np.linalg.matrix_rank(subs[[i for i in ok]]) if len(ok) else 0
    print("h=t^2+%s t+%s"%(a,b),"zero-sum sets:",len(good),"rank",rk)
