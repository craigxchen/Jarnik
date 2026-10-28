# Corollary check: 3-subtuples of the Thm 4.1 family, after their own gcd reduction, sit on arcs ~ R_3^{1/3}
import sys; sys.path.insert(0,'.')
from gi import *
exec(open('k4_lattice.py').read().split('cuts=')[0])
from itertools import combinations
from math import log, sqrt
for s in [10**6,10**12,10**20,10**30]:
    Z=pts(s); g=ggcdlist(Z); zeta=[gexact(z,g) for z in Z]
    for drop in range(4):
        T=[zeta[i] for i in range(4) if i!=drop]; g3=ggcdlist(T); eta=[gexact(z,g3) for z in T]
        R3=sqrt(gnorm(eta[0])); d=max(sqrt(gnorm((a[0]-b[0],a[1]-b[1]))) for a,b in combinations(eta,2))
        print(s,drop,'exp=%.4f'%(log(d)/log(R3)),'diam/R3^(1/3)=%.3f'%(d/R3**(1/3)))
