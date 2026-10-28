import sys, random
from math import log
sys.path.insert(0,'.')
import ref_phiN as P
from ref_lattice import lattice
import flint
x=P.x
Q1,Q2=P.Qj(0),P.Qj(1); D=P.add(Q2,P.neg(Q1))
for N in map(int,sys.argv[1:]):
    R=P.add(Q1,P.mul(N,D)); a,b,c=P.lifts(R)[0]
    d,ok,cop=P.profile([x[0],x[1],x[2],a,b,c]); w=sum(log(v) for v in d.values())/25/2
    val=lambda C: x[min(C)] if min(C)<=2 else x[3]
    opts=[]
    for C in P.sides:
        dd=d[C]
        if dd<=1: continue
        fac=[(int(p),e) for p,e in flint.fmpz(dd).factor()]
        divs=[1]
        for p,e in fac: divs=[q*p**k for q in divs for k in range(e+1)]
        L=log(dd); win=[q for q in divs if 0.25*L<=log(q)<=0.75*L]
        if win: opts.append((C,win))
    rnd=random.Random(12345); mx=-1e9; tot=0; ns=1000
    for it in range(ns):
        conds=[]
        for C,win in opts:
            e=rnd.choice(win); v=val(C); conds.append((v.numerator,v.denominator,e))
        M,vec=lattice(conds)
        sk=(log(M)/2-log(max(abs(vec[0]),abs(vec[1]))))/w
        mx=max(mx,sk); tot+=sk
    k=len(opts)
    print(N,'w=%.2f'%w,'contacts used',k,'max skew/w %.3f mean %.3f'%(mx,tot/ns),
          'needed at eta=0 (k/2-8): %.2f'%(k/2-8),'needed at eta=1/5 (k(1-eta)/2-8(1+eta)): %.2f'%(k*0.4-9.6),flush=True)
