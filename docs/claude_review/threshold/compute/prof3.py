import time, sys
from combi import *
import iso, modla
M,pp,nn=map(int,sys.argv[1:4]); lam=tuple(map(int,sys.argv[4].split(',')))
orbs=shell_orbits(M,pp,nn)
for dg in (None, phi_floor(M,pp,nn)):
    T=time.time()
    r=iso.lam_profile(lam,orbs,modla.P1,dguess=dg)
    print('dguess',dg,'N',r['N'],'reg',r['reg'],'ngen',r['ngen'],'time %.2f'%(time.time()-T), r['profile'][-4:])
