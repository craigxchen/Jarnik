from fractions import Fraction as F
import itertools
from collections import Counter
q=[F(-9,4),F(-2,3),F(1,4),F(8,3)]; eps=[1,-1,-1,1]
roots=[]
for qj,e in zip(q,eps):
    r=(e/qj, qj)  # x+iy
    for s in (1,-1):
        x,y=s*r[0],s*r[1]
        if y<0: x,y=x,-y  # conjugate into UHP (conj of -r = -rbar)
        roots.append((x,y))
print(roots); assert len(set(roots))==8 and all(y>0 for _,y in roots)
def mom(x,y): return (y, 2*x*y, 3*x*x*y-y**3)  # Im rho, Im rho^2, Im rho^3
H=[mom(*r) for r in roots]
vecs=[m for m in itertools.product(range(-2,3),repeat=8) if sum(map(abs,m))==4 and all(sum(mi*h[c] for mi,h in zip(m,H))==0 for c in range(3))]
tern=[m for m in vecs if max(map(abs,m))==1]
print('weight-4 kernel vectors (both signs):',len(vecs),' ternary:',len(tern))
# full weight<=4 check incl. weight 2,3
for w in (1,2,3):
    assert not [m for m in itertools.product(range(-3,4),repeat=8) if sum(map(abs,m))==w and all(sum(mi*h[c] for mi,h in zip(m,H))==0 for c in range(3))]
code=[15,51,60,85,170,195,204,240]
c=Counter(bin(a^b).count('1') for a,b in itertools.combinations(code,2)); print('code edge weights',c)
sup=Counter(a^b for a,b in itertools.combinations(code,2) if bin(a^b).count('1')==4); print(len(sup), set(sup.values()))
