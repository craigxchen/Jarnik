# Defect sizes on Astra's 191-digit four-star instance (unsheared columns), measured against P=N^(1/10).
from fractions import Fraction
from itertools import combinations
from math import prod, log
roots=(1,-9,-7,2,4); marks=[r*r for r in roots]; edges=list(combinations(range(5),2))
scale=6_350_400; parameter=3_175_200_023; t=Fraction(parameter,scale)
d_real=t**4+25*t**2+504; d_imag=10*t**3+250*t
cols=[]
for u in marks:
    mult=Fraction(u-1,u)*(u+t**2)
    re=scale**4*(d_real+mult*t/10); im=scale**4*(d_imag+mult)
    cols.append((int(re),int(im)))
blocks={e:parameter**2+(scale*roots[e[0]]*roots[e[1]])**2 for e in edges}
N=prod(blocks.values()); P=N**0.1 if False else None
lnN=log(N); lnP=lnN/10
for i in range(5):
    n=prod(b for e,b in blocks.items() if i in e); x,y=cols[i]; d=x*x+y*y-n
    iso=[e for e,b in blocks.items() if i in e and (x*x+y*y)%b==0]
    print('mark',i+1,'log|delta|/logP =', 'exact' if d==0 else round(log(abs(d))/lnP,3), ' log n/logP=',round(log(n)/lnP,3),' isotropic blocks at mark:',[(a+1,b+1) for a,b in iso])
# interpolation identity (IV.5 / Part III (8)): sum_i n_i V_i / prod_{j!=i} m_ij = 0 ?
def m(i,j): return cols[i][0]*cols[j][1]-cols[i][1]*cols[j][0]
S=[Fraction(0),Fraction(0)]
for i in range(5):
    n=prod(b for e,b in blocks.items() if i in e); den=prod(m(i,j) for j in range(5) if j!=i)
    S[0]+=Fraction(n*cols[i][0],den); S[1]+=Fraction(n*cols[i][1],den)
print('linear identity sum n_i V_i/prod m_ij == 0 ?', S[0]==0 and S[1]==0)
