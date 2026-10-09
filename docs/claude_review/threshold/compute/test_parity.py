import time
from combi import *
from iso import *
for (M,pp) in [(3,1),(3,3),(4,2),(4,3),(5,3),(5,4),(6,3),(6,4),(7,3),(7,6)]:
    orbs=shell_orbits(M,pp,pp)
    for kind,orb in (('shell',orbs),('slice',slice_orbits(M,pp,pp))):
        t1=t2=0
        for lam in partitions(M):
            a=time.time(); r1=lam_profile(lam,orb,P1,dguess=phi_floor(M,pp,pp)); t1+=time.time()-a
            a=time.time(); r2=lam_profile_parity(lam,orb,P1,dguess=phi_floor(M,pp,pp)); t2+=time.time()-a
            if r1['N']==0: continue
            assert r1['reg']==r2['reg'] and r1['profile']==r2['profile'], (M,pp,kind,lam,r1['profile'],r2['profile'])
        print(M,pp,kind,'ok  time plain %.2f parity %.2f'%(t1,t2),flush=True)
