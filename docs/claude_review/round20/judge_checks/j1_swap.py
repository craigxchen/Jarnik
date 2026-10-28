# Judge check J1 (evidence; exact integer lattice data, float angles).
# F3 = smoothed-radius, sharp-angle sector counts with badly approximable (quadratic-irrational) edges.
# Killer's c3b: the MINIMAL graft of a 40-prime fake violates F3 (mean |err| 2.30 vs genuine 0.141 at X=1.5e6).
# Claim J1: a NON-multiplicative "swap extension" restores F3. For every S-smooth ideal d (moved to psi_d), pick a distinct
# non-S-smooth ideal alpha_d with genuine angle within eta of psi_d, in the same cell of the edge partition, and
# |N alpha_d - N d| minimal; give alpha_d the genuine angle theta(d). The multiset {(norm, angle)} then changes only by
# norm shifts, so each sector count changes by <= sup|f'| * sum |N alpha_d - N d| / X.
# Same fake (seed 7, s=40 primes in [1009,1601]) and same 30 sectors as the killer's c3b script.
import numpy as np, math, sympy, sys, os
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
slopes=sorted(set(math.sqrt(m)%1 for m in range(2,400) if math.isqrt(m)**2!=m))
uu=np.linspace(0.5,2,200001); fp=np.abs(np.gradient(f(uu),uu)).max()
print(f"S range {S[0]}..{S[-1]} s={len(S)}  sup|f'|={fp:.4f}",flush=True)
for X in [1_500_000, 2_500_000]:
    x,y,n,th,wt=run(X)
    spart=np.ones(len(x)); shift=np.zeros(len(x))
    for j,q in enumerate(S):
        aa,bb=gp(q); u=vals(x,y,q,aa,bb); w=vals(x,y,q,aa,-bb)
        spart*=float(q)**(u+w); shift+=(u-w)*(phi[j]-th_S[j])
    sm=np.isclose(spart,n.astype(float)); thmin=np.where(sm,np.mod(th+shift,math.pi/2),th)
    secs=[]
    for t in range(30):
        s1,s2=sorted(rng.choice(slopes,2,replace=False)); secs.append((math.atan(s1),math.atan(s2)))
    E=np.array(sorted(set(a for s in secs for a in s)))
    cell=np.searchsorted(E,th,side='left')  # cell index w.r.t. all edges (edges irrational: no ties)
    eta=X**-0.4
    order=np.argsort(th); ths=th[order]
    used=np.zeros(len(x),bool); used[sm]=True   # S-smooth ideals are never used as partners
    thsw=thmin.copy(); dn=[]; dth=[]; fail=0
    for i in np.nonzero(sm)[0]:
        psi=thmin[i]; c=np.searchsorted(E,psi,side='left')
        lo=np.searchsorted(ths,psi-eta); hi=np.searchsorted(ths,psi+eta)
        cand=order[lo:hi]; cand=cand[(~used[cand])&(cell[cand]==c)]
        if len(cand)==0: fail+=1; continue
        k=cand[np.argmin(np.abs(n[cand]-n[i]))]
        used[k]=True; thsw[k]=th[i]; dn.append(abs(int(n[k])-int(n[i]))); dth.append(abs(th[k]-psi))
    dn=np.array(dn); bound=fp*dn.sum()/X
    ge=[];gr=[];gs=[]
    for (a1,a2) in secs:
        main=(a2-a1)/2*X*If
        ge.append(float((wt*((th>=a1)&(th<=a2))).sum())-main)
        gr.append(float((wt*((thmin>=a1)&(thmin<=a2))).sum())-main)
        gs.append(float((wt*((thsw>=a1)&(thsw<=a2))).sum())-main)
    ge=np.abs(ge); gr=np.abs(gr); gs=np.abs(gs); dif=np.abs(np.array(gs)-np.array(ge))
    print(f"X={X} #S-smooth={int(sm.sum())} partner failures={fail} eta={eta:.2e} max|dn|={dn.max()} mean|dn|={dn.mean():.1f} "
          f"(X^0.5={X**0.5:.0f}) max|dtheta|={max(dth):.2e}", flush=True)
    print(f"   genuine |err| max={ge.max():.3f} mean={ge.mean():.3f} | minimal graft max={gr.max():.3f} mean={gr.mean():.3f} | "
          f"swap extension max={gs.max():.3f} mean={gs.mean():.3f} | max|swap-genuine|={dif.max():.4f} <= proved bound {bound:.4f}", flush=True)
