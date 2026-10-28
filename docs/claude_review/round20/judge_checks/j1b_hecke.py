# Judge check J1b (evidence): smooth all-ideal Hecke sums F_k(X)=sum f(N a/X) e^{4ik theta(a)}, k=1,2,3, for the genuine
# ensemble, the minimal graft, and a swap extension with eta = X^{-1/4} (partner chosen to minimise |dn|/X + |dtheta|).
import numpy as np, math, sympy, os
os.chdir('/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round20/kill_composite_checks')
exec(open('c3_ba_sector.py').read().split("print(f\"int f")[0])
rng=np.random.default_rng(7)
S=[q for q in sympy.primerange(1000,2400) if q%4==1][:40]
def gp(q):
    for bb in range(1,math.isqrt(q)+1):
        aa=math.isqrt(q-bb*bb)
        if aa*aa==q-bb*bb and aa>bb: return aa,bb
th_S=np.array([math.atan2(gp(q)[1],gp(q)[0]) for q in S]); phi=rng.uniform(0,math.pi/2,len(S))
def vals(x,y,q,aa,bb):
    v=np.zeros(len(x),dtype=np.int64); act=np.ones(len(x),bool); x=x.copy(); y=y.copy()
    while act.any():
        re_=x*aa+y*bb; im=y*aa-x*bb; ok=act&(re_%q==0)&(im%q==0)
        if not ok.any(): break
        x=np.where(ok,re_//q,x); y=np.where(ok,im//q,y); v+=ok; act=ok
    return v
for X in [1_500_000]:
    x,y,n,th,wt=run(X)
    spart=np.ones(len(x)); shift=np.zeros(len(x))
    for j,q in enumerate(S):
        aa,bb=gp(q); u=vals(x,y,q,aa,bb); w=vals(x,y,q,aa,-bb)
        spart*=float(q)**(u+w); shift+=(u-w)*(phi[j]-th_S[j])
    sm=np.isclose(spart,n.astype(float)); thmin=np.where(sm,np.mod(th+shift,math.pi/2),th)
    order=np.argsort(th); ths=th[order]; used=sm.copy(); thsw=thmin.copy()
    for eta in [X**-0.25, X**-0.4]:
        used=sm.copy(); thsw=thmin.copy()
        for i in np.nonzero(sm)[0]:
            psi=thmin[i]; lo=np.searchsorted(ths,psi-eta); hi=np.searchsorted(ths,psi+eta)
            cand=order[lo:hi]; cand=cand[~used[cand]]
            k=cand[np.argmin(np.abs(n[cand]-n[i])/X+np.abs(th[cand]-psi))]
            used[k]=True; thsw[k]=th[i]
        for kk in [1,2,3]:
            Fg=np.sum(wt*np.exp(4j*kk*th)); Fm=np.sum(wt*np.exp(4j*kk*thmin)); Fs=np.sum(wt*np.exp(4j*kk*thsw))
            print(f"X={X} eta={eta:.2e} k={kk}: |F| genuine={abs(Fg):.3e} minimal graft={abs(Fm):.3e} swap={abs(Fs):.3e}",flush=True)
