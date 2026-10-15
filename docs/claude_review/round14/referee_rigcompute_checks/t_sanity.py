from reflib import *
ps=[5,13,17,29,37,53,97,353]; pis=[gprime(p) for p in ps]
bad=0; cnt=0
for q in primerange(3,200):
    if q in ps: continue
    L,iota,m=genuine_L(ps,q)
    assert [x%2 for x in L]==nu_par(ps,q),(q,L)
    for v in monos(len(ps),2):
        d=sum(a*b for a,b in zip(v,L))%m
        z=Zv(pis,v); D=Dvals(z)
        for s in range(4):
            col=(d==(s*iota)%m); dv=(D[s]%q==0)
            cnt+=1
            if col!=dv: bad+=1
print("parity OK; collision<->D_s mismatches",bad,"of",cnt)
