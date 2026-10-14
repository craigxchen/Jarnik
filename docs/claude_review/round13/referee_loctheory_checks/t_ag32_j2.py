# AG(3,2) image at j=2: does J(P) leave the grid below degree 4j=8?  full column rank mod p of the reduced
# evaluation matrix in degree <= e certifies h_P(e) = N_{<=e} exactly (hence J(P)_{<=e} = grid).
import itertools, sys, flint, time
from rlib import *
j=2; b=8; t=4
V=[a for a in itertools.product((0,1),repeat=t) if any(a)]
F=[a for a in V if sum(a)%2==1]
codes=[s for s in itertools.product((0,1),repeat=b) if any(s) and all(sum(F[i][q] for i in range(b) if s[i])%2==0 for q in range(t))]
P=[y for y in itertools.product(range(-j,j+1),repeat=b) if sum(y)%4==0 and all(sum(y[i] for i in range(b) if s[i])%2==0 for s in codes)]
print("|P| =",len(P)); sys.stdout.flush()
cap=[5]*b
for e in range(6,10):
    T0=time.time()
    N=len(monos_le(b,e,cap)); r=h_func_modp(P,e,cap)
    print("e",e,"#reduced monos<=e",N,"mod-p rank",r,"full" if r==N else "NOT full (evidence of non-grid J in degree <= e)","%.0fs"%(time.time()-T0)); sys.stdout.flush()
    if r<N: break
