import numpy as np, time, sys
from combi import *
from iso import *
from direct import direct_profile
p=P1
# 1. method C vs method R agreement on random points
rng=np.random.default_rng(5)
for M in range(2,7):
    for lam in partitions(M):
        ld=LamData(lam)
        K=rng.integers(-5,6,(50,M))
        m0=ld.method
        ld.method='C'; ld.perms,ld.signs=group_perms(ld.cols,M,True); FC=ld.F_values(K,p)
        ld.method='R'; ld.perms,_=group_perms(ld.rows,M,False)
        inv=np.empty_like(ld.perms)
        for i,pr in enumerate(ld.perms): inv[i,pr]=np.arange(M)
        ld.perms_inv=inv; FR=ld.F_values(K,p)
        assert np.array_equal(FC,FR),lam
print('method C == method R ok')
# 2. totals vs direct Hilbert function
def check(M,pp,nn):
    orbs=shell_orbits(M,pp,nn)
    pts=shell_points(M,pp,nn)
    dp=direct_profile(pts,p)
    tot=None; regs={}
    for lam in partitions(M):
        r=lam_profile(lam,orbs,p)
        if r['N']==0: continue
        regs[lam]=r['reg']
        prof=r['profile']
        L=max(len(dp),len(prof))
        pr=[prof[min(d,len(prof)-1)] for d in range(L)]
        contrib=[hook_dim(lam)*x for x in pr]
        if tot is None: tot=contrib
        else:
            L2=max(len(tot),L); tot=[ (tot[min(d,len(tot)-1)] if True else 0)+contrib[min(d,len(contrib)-1)] for d in range(L2)]
    dpp=[dp[min(d,len(dp)-1)] for d in range(len(tot))]
    ok = tot==dpp
    print(M,pp,nn,'|S|=',len(pts),'direct reg',len(dp)-1,'iso max reg',max(regs.values()),'profiles equal',ok, flush=True)
    if not ok: print(' direct',dp,'\n iso   ',tot)
    return ok
for M in range(2,6):
    for t in range(1,9):
        for s in range(t%2,t+1,2):
            pp=(t+s)//2; nn=(t-s)//2
            check(M,pp,nn)
