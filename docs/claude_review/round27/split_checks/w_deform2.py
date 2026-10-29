# Continuation away from the W-curve: prescribe s = (sm+e, sm-e, -sm+i e, -sm-i e) and track the curve as e grows;
# report the dominance singular values (rank 4 <=> the 4-branch family dominates T).
import numpy as np, itertools, sys
exec(open('w_deform.py').read().split("# deform: split the cusps")[0])
def eqs_fixed(w,si):
    tp,lam=w[:nP],np.concatenate([w[nP:nP+4],[1.0]]); out=[]
    for x in si:
        u=U(tp,lam,x); out.extend(u[:4]/u[4]-1)
    return np.array(out)
def jac_fixed(w,si,h=1e-8):
    f0=eqs_fixed(w,si); J=np.zeros((len(f0),len(w)),dtype=complex)
    for k in range(len(w)):
        e=np.zeros(len(w),dtype=complex); e[k]=h; J[:,k]=(eqs_fixed(w+e,si)-f0)/h
    return J
w=np.concatenate([np.array(places),lam0[:4]])
for e in [0.0,0.01,0.02,0.05,0.1,0.2,0.3,0.5,0.8]:
    si=np.array([sm+e,sm-e,-sm+1j*e,-sm-1j*e])
    with np.errstate(all='ignore'):
        for it in range(100):
            f=eqs_fixed(w,si)
            if np.linalg.norm(f)<1e-13: break
            w=w+np.linalg.lstsq(jac_fixed(w,si),-f,rcond=None)[0]
    res=np.linalg.norm(eqs_fixed(w,si))
    wf=np.concatenate([w,si])
    J=jac(wf); sv=np.linalg.svd(J,compute_uv=False); rk=int(np.sum(sv>1e-7*sv[0]))
    _,_,Vh=np.linalg.svd(J); K=Vh[rk:].conj().T
    x0=0.3+0.2j
    y=lambda ww,x: (lambda tp,lam,_: U(tp,lam,x)[:4]/U(tp,lam,x)[4])(*unpack(ww))
    h=1e-6; cols=[(y(wf+h*K[:,k],x0)-y(wf-h*K[:,k],x0))/(2*h) for k in range(K.shape[1])]+[(y(wf,x0+h)-y(wf,x0-h))/(2*h)]
    dsv=np.linalg.svd(np.array(cols).T,compute_uv=False)
    print(f"e={e}: residual {res:.1e}, solution-set dim {len(wf)-rk}, dominance sv/sv0 = {np.array2string(dsv/dsv[0],precision=3)}",flush=True)
