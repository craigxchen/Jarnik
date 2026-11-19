from pell_spec import *
from fractions import Fraction as F
E5=[(0,1),(0,2),(1,2),(0,3)]
for n in [1,2]:
    for T in [3,7]:
        seeds=[F(1),F(3),F(5),F(8)]
        sv=jac_rank(n,seeds,E5,T)
        print("n",n,"T",T,"singular values of d(args)/d(seeds):",sv)
# a degenerate check: drop edge 14 -> 4 marks with r=4 should give rank 3
E4=[(0,1),(0,2),(1,2)]
print("4 marks rank check:",jac_rank(1,[F(1),F(3),F(5),F(8)],E4,3))
