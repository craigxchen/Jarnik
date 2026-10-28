# F3 check: a position-dependent composite input with precision X^{o(1)}.
# S_I(X) = sum_{ideals alpha} f(N alpha / X) 1[vartheta(alpha) in I],  f = C^infty bump on (1/2,2), I=[psi1,psi2].
# Main term (|I|/2) X int f.  Error for BA edges (quadratic-irrational slopes) vs rational edges.
# Then the change caused by the MINIMAL graft of a fake with S = 12 primes near 10^3 (Haar-random fake angles).
import numpy as np, math
from scipy.integrate import quad
def f(u):
    u=np.asarray(u,float); out=np.zeros_like(u); m=(u>0.5)&(u<2)
    out[m]=np.exp(-1.0/((u[m]-0.5)*(2-u[m])))
    return out
If=quad(lambda u: float(f(np.array([u]))[0]),0.5,2,epsabs=1e-14,epsrel=1e-13)[0]
edges={'BA (sqrt2-1, (sqrt5-1)/2)':(math.atan(math.sqrt(2)-1),math.atan((math.sqrt(5)-1)/2)),
       'BA (sqrt7-2, sqrt3-1)':(math.atan(math.sqrt(7)-2),math.atan(math.sqrt(3)-1)),
       'rational (1/3, 1/2)':(math.atan(1/3),math.atan(1/2)),
       'rational (2/5, 3/5)':(math.atan(2/5),math.atan(3/5))}
def gen(X):
    r=math.isqrt(2*X); 
    for x0 in range(1,r+1,2000):
        xs=np.arange(x0,min(x0+2000,r+1),dtype=np.int64)
        for x in xs:
            ylo=math.isqrt(max(X//2-int(x)*int(x),0)); yhi=math.isqrt(max(2*X-int(x)*int(x),0))
            if 2*X-int(x)*int(x)<0: continue
            y=np.arange(max(ylo-1,0),yhi+1,dtype=np.int64)
            yield int(x),y
def run(X):
    xs=[];ys=[]
    for x,y in gen(X):
        xs.append(np.full(len(y),x,dtype=np.int64)); ys.append(y)
    x=np.concatenate(xs); y=np.concatenate(ys); n=x*x+y*y
    keep=(n>X/2)&(n<2*X); x=x[keep]; y=y[keep]; n=n[keep]
    th=np.arctan2(y.astype(float),x.astype(float)); wt=f(n/X)
    return x,y,n,th,wt
print(f"int f = {If:.12f}")
for X in [10**4,10**5,10**6,4*10**6,10**7]:
    x,y,n,th,wt=run(X)
    line=f"X={X:9d}"
    for name,(a1,a2) in edges.items():
        s=float((wt*((th>=a1)&(th<=a2))).sum()); main=(a2-a1)/2*X*If
        line+=f" | {name}: err={s-main:+9.3f}"
    print(line,flush=True)

# ---- minimal graft of a fake with S = 12 primes near 1000, Haar-random fake angles (seeded) ----
import sympy
rng=np.random.default_rng(20)
S=[q for q in sympy.primerange(1000,1400) if q%4==1][:12]
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
print("S =",S)
for X in [10**6,4*10**6,10**7]:
    x,y,n,th,wt=run(X)
    spart=np.ones(len(x)); shift=np.zeros(len(x))
    for j,q in enumerate(S):
        aa,bb=gp(q); u=vals(x,y,q,aa,bb); w=vals(x,y,q,aa,-bb)
        spart*=float(q)**(u+w); shift+=(u-w)*(phi[j]-th_S[j])
    sm=np.isclose(spart,n.astype(float))
    thmin=np.where(sm,np.mod(th+shift,math.pi/2),th)
    line=f"X={X:9d} #S-smooth in window={int(sm.sum())}"
    for name,(a1,a2) in edges.items():
        if not name.startswith('BA'): continue
        g=float((wt*((th>=a1)&(th<=a2))).sum()); m=float((wt*((thmin>=a1)&(thmin<=a2))).sum()); main=(a2-a1)/2*X*If
        line+=f" | {name}: genuine err={g-main:+.3f} minimal-graft err={m-main:+.3f}"
    print(line,flush=True)
