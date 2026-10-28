import sys, itertools; sys.path.insert(0,sys.argv[1])
from regp import T
from rec_check_funcs import F
from fractions import Fraction
ok=True
for a in range(6):
  for b in range(6):
    S=T(5,a,b)
    for A in itertools.combinations(range(5),2):
      v=[sum(s[j] for j in A) for s in S]
      if (max(v),min(v))!=(a,-b): ok=False
print("Cor 3.1: every pair width of T_5(a,b) = a+b with range [-b,a] (a,b<=5):",ok)
eq=[(a,b) for a in range(200) for b in range(200) if 3*F(5,a,b)==5*(a+b)]
print("F5<=5(a+b)/3 on [0,200)^2:",all(3*F(5,a,b)<=5*(a+b) for a in range(200) for b in range(200)),
      "; equality set == {(3s,3s)}:", eq==[(3*s,3*s) for s in range(67)])
