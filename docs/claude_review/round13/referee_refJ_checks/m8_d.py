from rlib import *
import pickle, random
A,H,dimJ,W=pickle.load(open('m8_a.pkl','rb'))
fam,done=pickle.load(open('m8_c.pkl','rb'))
p=P; random.seed(3)
def fval_modp(f,V,pt):
    ys=[sum(a*b for a,b in zip(v,pt))%p for v in V]
    s=0
    for ex,c in f.items():
        t=c%p
        for e_,y in zip(ex,ys):
            if e_: t=t*pow(y,e_,p)%p
        s=(s+t)%p
    return s
def rows_for(deg,pts,which):
    rows=[]
    pm={}
    def mono_vals(e):
        if e not in pm: pm[e]=[[evmono(m,pt)%p for pt in pts] for m in monos(7,e)]
        return pm[e]
    for v,t,h in W:
        if h+1<=deg:
            fv=[pow(sum(a*b for a,b in zip(v,pt))%p,h+1,p) for pt in pts]
            for mv in mono_vals(deg-h-1): rows.append([x*y%p for x,y in zip(fv,mv)])
    for key in fam:
        if key[0] not in which: continue
        V=fam[key]; g=done[key]
        # minimal set: generators in degree e not in R_1*J_{e-1} not computed; use all degrees e (redundant but valid)
        for e in range(deg+1):
            for f in g[e]:
                fv=[fval_modp(f,V,pt) for pt in pts]
                for mv in mono_vals(deg-e): rows.append([x*y%p for x,y in zip(fv,mv)])
    return rows
for which in [(), ('pair',), ('nonline',), ('line',), ('pair','line','nonline')]:
    out=[]
    for deg in (5,6,7,8):
        N=math.comb(deg+6,6)
        pts=[tuple(random.randrange(p) for _ in range(7)) for _ in range(N+10)]
        rows=rows_for(deg,pts,which)
        out.append(nmod_mat(rows,p).rank())
    print("J_w +",which,"mod-p eval ranks deg5..8:",out,"  dim J(A):",dimJ[5:9],flush=True)
