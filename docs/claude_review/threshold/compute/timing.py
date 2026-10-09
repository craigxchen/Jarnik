import sys, time
from combi import *
from iso import *
M,pp,nn=map(int,sys.argv[1:4])
orbs=shell_orbits(M,pp,nn)
T=time.time()
for lam in partitions(M):
    r=lam_profile(lam,orbs,P1)
    print(lam,'N',r['N'],'reg',r['reg'],'meth',r['method'],'time %.2f'%r['time'], 'tF %.2f'%r.get('t_F',0), flush=True)
print('total',time.time()-T)
