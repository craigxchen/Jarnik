"""(a) rank of the block-argument -> anchored-argument map, and the Haar arc fraction."""
import itertools, random, math
from fractions import Fraction
import sympy as sp
pairs=list(itertools.combinations(range(5),2))
# arg z_k = sum_e theta_e - 2 sum_{e ni k} theta_e ; anchored: arg z_k - arg z_0
M=sp.Matrix([[ -2*((k in e)-(0 in e)) for e in pairs] for k in range(1,5)])
print("rank of 4x10 anchored map:",M.rank())
# index of the image lattice in Z^4 (Smith form): image of Z^10 under M
from sympy.matrices.normalforms import smith_normal_form
print("Smith form diag:",[smith_normal_form(M,domain=sp.ZZ)[i,i] for i in range(4)])
# Haar fraction: 5 iid uniform marks on circle all in an arc of angular length d
random.seed(1)
def inarc(th,d):
    th=sorted(th); gaps=[th[i+1]-th[i] for i in range(4)]+[th[0]+2*math.pi-th[4]]
    return 2*math.pi-max(gaps)<=d
for d in [0.3,1.0,2.0,3.0]:
    n=400000; c=sum(inarc([random.uniform(0,2*math.pi) for _ in range(5)],d) for _ in range(n))
    print("delta=%.1f MC %.5f formula 5(d/2pi)^4 = %.5f"%(d,c/n,5*(d/2/math.pi)**4))
# pushforward through the actual map from ten uniform block arguments
for d in [1.0,2.0]:
    n=200000;c=0
    for _ in range(n):
        th=[random.uniform(0,2*math.pi) for _ in range(10)]
        ang=[0.0]+[sum(int(M[k,j])*th[j] for j in range(10))%(2*math.pi) for k in range(4)]
        c+=inarc(ang,d)
    print("via block map delta=%.1f MC %.5f formula %.5f"%(d,c/n,5*(d/2/math.pi)**4))
