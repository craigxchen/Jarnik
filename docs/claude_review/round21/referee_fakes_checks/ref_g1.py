import math, json, numpy as np
from ref_seq import enum
HP=math.pi/2
def g1(p,nmax=60):
    # forbidden open arcs for n t: centers k*HP/n (axis) and (pi/4 + k HP)/n (diag), radius r/n
    bad=[]
    for n in range(1,nmax+1):
        ra=math.asin(p**(-n/2)); rd=math.asin(p**(-n/2)/math.sqrt(2))
        if ra/n<1e-15: break
        for k in range(n):
            c=k*HP/n; bad.append((c-ra/n,c+ra/n))
            c=(math.pi/4+k*HP)/n; bad.append((c-rd/n,c+rd/n))
    ivs=[]
    for lo,hi in bad:
        ivs+= [(lo,hi),(lo+HP,hi+HP),(lo-HP,hi-HP)]
    ivs.sort(); covered=0.0; cur_lo,cur_hi=None,None; comps=0
    merged=[]
    for lo,hi in ivs:
        if merged and lo<=merged[-1][1]: merged[-1][1]=max(merged[-1][1],hi)
        else: merged.append([lo,hi])
    tot=0.0
    for lo,hi in merged:
        a=max(lo,0.0); b=min(hi,HP)
        if b>a: tot+=b-a
    gaps=sum(1 for i in range(len(merged)-1) if merged[i][1]<HP and merged[i+1][0]>0 and merged[i][1]>=0)
    return 1-tot/HP, gaps
for p in [13,17]:
    m,c=g1(p); print(json.dumps(dict(p=p,mu1=round(m,5),components_approx=c)))
# table columns
rows=[[13],[13,17],[13,17,29],[13,17,29,37],[13,17,29,37,41],[97,101,109,113,137,149],[5]]
for ps in rows:
    print(json.dumps(dict(primes=ps,Lambda=round(sum(p**-0.5 for p in ps),3),Pi=round(float(np.prod([1+2/(math.sqrt(p)-1) for p in ps])),3))))
# MC for {13,17}, T=26
rng=np.random.default_rng(12345)
ps=[13,17]; logs=[math.log(p) for p in ps]; V=enum(logs,26); w=np.abs(V)@np.array(logs)
ta=np.arcsin(np.exp(-w/2)); td=np.arcsin(np.exp(-w/2)/math.sqrt(2))
n=100000; X=rng.random((n,2))*HP; ok=np.ones(n,bool)
for s in range(0,n,5000):
    t=X[s:s+5000]@V.T
    ra=np.mod(t,HP); da=np.minimum(ra,HP-ra); rd=np.mod(t-math.pi/4,HP); dd=np.minimum(rd,HP-rd)
    ok[s:s+5000]=((da>=ta)&(dd>=td)).all(1)
print(json.dumps(dict(mc="13,17",T=26,nV=len(V),frac=float(ok.mean()))))
