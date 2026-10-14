from rlib import *
import pickle, time, random
A,H,dimJ,W=pickle.load(open('m8_a.pkl','rb'))
random.seed(1)
p=P
def evalpoly_modp(f,pt):
    s=0
    for ex,c in f.items():
        t=c
        for e_,y in zip(ex,pt):
            if e_: t=t*pow(y,e_,p)%p
        s=(s+t)%p
    return s
# generators of J_w as (degree, eval function). width gen: l_v^{h+1}
def gen_rows_modp(gens, deg, pts):
    """gens: list of (e, callable value at point). rows: g*m for monomial m of degree deg-e, evaluated at pts"""
    rows=[]
    for e,fval in gens:
        if e>deg: continue
        fv=[fval(pt) for pt in pts]
        for m in monos(7,deg-e):
            rows.append([fv[i]*evmono(m,pts[i])%p for i in range(len(pts))])
    return rows
def lin_val(v): return lambda pt: pow(sum(a*b for a,b in zip(v,pt))%p, 1, p)
wgens=[(h+1,(lambda v,h: (lambda pt: pow(sum(a*b for a,b in zip(v,pt))%p,h+1,p)))(v,h)) for v,t,h in W]
res={}
for deg in (6,7,8):
    N=math.comb(deg+6,6)
    pts=[tuple(random.randrange(p) for _ in range(7)) for _ in range(N+20)]
    rows=gen_rows_modp(wgens,deg,pts)
    r=nmod_mat(rows,p).rank()
    res[deg]=r
    print("deg",deg,"#rows",len(rows),"modp eval-rank J_w (lower bound on Q dim)",r,"dim J(A)",dimJ[deg],flush=True)
pickle.dump(res,open('m8_b.pkl','wb'))
# exact over Q in degree 6 and 7: coefficient matrix
for deg in (6,7):
    basis=monos(7,deg); idx={m:i for i,m in enumerate(basis)}
    rows=[]
    for v,t,h in W:
        if h+1>deg: continue
        f=ppow(linform(v),h+1)
        for m in monos(7,deg-h-1):
            rows.append(vec(pmul(f,{m:1}),idx))
    t0=time.time(); r=fmpz_mat(rows).rank()
    print("deg",deg,"EXACT Q rank of (J_w)_deg:",r,"time",round(time.time()-t0,1),flush=True)
