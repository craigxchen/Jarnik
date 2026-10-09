from combi import *
import sys
c=1.05e-9
def Ns(M,p,n):
    O=shell_orbits(M,p,n); mus=[mu_of(v) for v in O]
    return [sum(kostka(l,m) for m in mus) for l in partitions(M)]
for M in range(3,10):
    cum=0
    for t in range(1,25):
        tt=0
        for s in range(t%2,t+1,2):
            p=(t+s)//2;n=(t-s)//2
            tt+=sum(c*x**3+0.05 for x in Ns(M,p,n))
        cum+=tt
        print(M,t,'shell-t time %.0f s, cumulative %.0f s'%(tt,cum),flush=True)
        if cum>40000: break
