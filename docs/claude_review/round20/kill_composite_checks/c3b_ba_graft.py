# F3 (evidence): over 30 sectors with quadratic-irrational edge slopes, compare the genuine error of the smoothed
# sector count with the change made by the minimal graft of a fake with s = 40 primes in [1000, 2400], X ~ P0^2.
import numpy as np, math, sympy
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
# quadratic-irrational slopes in (0,1): frac(sqrt(m)) for non-square m
slopes=sorted(set(math.sqrt(m)%1 for m in range(2,400) if math.isqrt(m)**2!=m))
print("S range",S[0],S[-1],"s=",len(S))
for X in [1_500_000, 2_500_000]:
    x,y,n,th,wt=run(X)
    spart=np.ones(len(x)); shift=np.zeros(len(x))
    for j,q in enumerate(S):
        aa,bb=gp(q); u=vals(x,y,q,aa,bb); w=vals(x,y,q,aa,-bb)
        spart*=float(q)**(u+w); shift+=(u-w)*(phi[j]-th_S[j])
    sm=np.isclose(spart,n.astype(float)); thmin=np.where(sm,np.mod(th+shift,math.pi/2),th)
    ge=[];gr=[]
    for t in range(30):
        s1,s2=sorted(rng.choice(slopes,2,replace=False)); a1,a2=math.atan(s1),math.atan(s2)
        main=(a2-a1)/2*X*If
        ge.append(float((wt*((th>=a1)&(th<=a2))).sum())-main); gr.append(float((wt*((thmin>=a1)&(thmin<=a2))).sum())-main)
    ge=np.abs(ge); gr=np.abs(gr)
    print(f"X={X} #S-smooth in window={int(sm.sum())} weight={float(wt[sm].sum()):.2f} | genuine |err|: max={ge.max():.3f} mean={ge.mean():.3f} | minimal graft |err|: max={gr.max():.3f} mean={gr.mean():.3f}",flush=True)
