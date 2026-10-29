# Dominance test by hitting random targets: from a 4-branch curve near W, solve for (curve in the 4-point family, s)
# with y(s) = y* for random targets y* near the curve.  Residual -> 0 for random y*  <=>  family dominates T.
# Control: the same with the cusp constraint (W-curves only): targets off W should NOT be hit.
import numpy as np, itertools, sys
exec(open('w_deform.py').read().split("# deform: split the cusps")[0])
rng2=np.random.default_rng(7)
def full_eqs(w,ystar,cusp=False):
    tp,lam=w[:nP],np.concatenate([w[nP:nP+4],[1.0]]); si=w[nP+4:nP+8]; s=w[nP+8]
    out=[]
    pts=si if not cusp else np.array([si[0],si[0]+1e-30,si[2],si[2]])
    for x in si:
        u=U(tp,lam,x); out.extend(u[:4]/u[4]-1)
    if cusp:  # derivative conditions: d/dx (u_j/u_5) = 0 at si[0], si[2]; and si[1]=si[0], si[3]=si[2]
        h=1e-6
        for x in (si[0],si[2]):
            dy=(U(tp,lam,x+h)[:4]/U(tp,lam,x+h)[4]-U(tp,lam,x-h)[:4]/U(tp,lam,x-h)[4])/(2*h); out.extend(dy)
        out.extend([si[1]-si[0],si[3]-si[2]])
    u=U(tp,lam,s); out.extend(u[:4]/u[4]-ystar)
    return np.array(out)
def jacf(w,ystar,cusp,h=1e-8):
    f0=full_eqs(w,ystar,cusp); J=np.zeros((len(f0),len(w)),dtype=complex)
    for k in range(len(w)):
        e=np.zeros(len(w),dtype=complex); e[k]=h; J[:,k]=(full_eqs(w+e,ystar,cusp)-f0)/h
    return J
# starting curve: the e=0.3 four-branch curve
def eqs_fixed(w,si):
    tp,lam=w[:nP],np.concatenate([w[nP:nP+4],[1.0]]); out=[]
    for x in si:
        u=U(tp,lam,x); out.extend(u[:4]/u[4]-1)
    return np.array(out)
w=np.concatenate([np.array(places),lam0[:4]]); e=0.3
si=np.array([sm+e,sm-e,-sm+1j*e,-sm-1j*e])
for it in range(100):
    f=eqs_fixed(w,si)
    if np.linalg.norm(f)<1e-13: break
    J=np.array([(eqs_fixed(w+1e-8*np.eye(len(w))[k],si)-f)/1e-8 for k in range(len(w))]).T
    w=w+np.linalg.lstsq(J,-f,rcond=None)[0]
for cusp in (False,True):
    wc0=np.concatenate([w,si,[0.3+0.2j]]) if not cusp else np.concatenate([np.array(places),lam0[:4],[sm,sm,-sm,-sm],[0.3+0.2j]])
    y0=(lambda ww: U(ww[:nP],np.concatenate([ww[nP:nP+4],[1.0]]),ww[nP+8])[:4]/U(ww[:nP],np.concatenate([ww[nP:nP+4],[1.0]]),ww[nP+8])[4])(wc0)
    hits=0; resids=[]
    for trial in range(8):
        ystar=y0+0.02*(rng2.normal(size=4)+1j*rng2.normal(size=4))
        wc=wc0.copy()
        with np.errstate(all='ignore'):
            for it in range(200):
                f=full_eqs(wc,ystar,cusp)
                if np.linalg.norm(f)<1e-12 or not np.all(np.isfinite(f)): break
                wc=wc+np.linalg.lstsq(jacf(wc,ystar,cusp),-f,rcond=None)[0]
        r=np.linalg.norm(full_eqs(wc,ystar,cusp)); resids.append(r); hits+= r<1e-9
    print("cusp-constrained (W-curves)" if cusp else "four-branch family","targets hit:",hits,"/8; residuals",np.array2string(np.array(resids),precision=1),flush=True)
