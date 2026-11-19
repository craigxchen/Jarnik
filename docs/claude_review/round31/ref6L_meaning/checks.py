# ref6L_meaning: independent small checks
import itertools, sympy as sp
from fractions import Fraction
# 1. VI.6 / VI.4 cut algebra: 5 marks, 15 cuts. d_ij = sum of weights of cuts separating i,j.
marks=range(5)
cuts=[frozenset([i]) for i in marks]+[frozenset(p) for p in itertools.combinations(marks,2)]
w=sp.symbols('w0:15')
d={}
for i,j in itertools.combinations(marks,2):
    d[(i,j)]=sum(w[c] for c,A in enumerate(cuts) if (i in A)!=(j in A))
W=sum(w)
Wsing=sum(w[:5])
print("sum d_ij == 6W-2Wsing:", sp.simplify(sum(d.values())-(6*W-2*Wsing))==0)
# if all d_ij = 3/5 W and Wsing=0: solve for pair weights
sol=sp.solve([d[k]-sp.Rational(3,5)*W for k in d]+[w[c] for c in range(5)]+[W-1],w,dict=True)
print("solution with d_ij=3W/5, singletons 0:",sol)
# LP sanity: min over nonneg weights (W=1) of max... check min_{w>=0, W=1} max_ij(-d_ij) i.e. max min d_ij = 3/5
import scipy.optimize as so
import numpy as np
A=[]; 
for k in d:
    row=[1.0 if ((k[0] in A_)!=(k[1] in A_)) else 0.0 for A_ in cuts]
    A.append(row)
A=np.array(A)
# maximize t s.t. A w >= t, sum w =1, w>=0
c=np.zeros(16); c[-1]=-1
Aub=np.hstack([-A,np.ones((10,1))]); bub=np.zeros(10)
Aeq=np.array([[1.0]*15+[0.0]]); beq=[1.0]
r=so.linprog(c,A_ub=Aub,b_ub=bub,A_eq=Aeq,b_eq=beq,bounds=[(0,None)]*15+[(None,None)])
print("max min_ij d_ij / W =",-r.fun, " -> width exponent", -r.fun/2, "(Liouville/CC: 3/10)")
# 2. Astra's v_n relation: exact, no proper vanishing subsum
n=sp.Symbol('n')
for N in [10,100,1000]:
    v=lambda m: (m+sp.sqrt(m*m-4))/2
    terms=[v(N+1),1/v(N+1),-v(N),-1/v(N),-1]
    tot=sp.nsimplify(sp.simplify(sum(terms)))
    zs=[]
    for r_ in range(1,5):
        for S in itertools.combinations(range(5),r_):
            val=sp.N(sum(terms[i] for i in S),50)
            if abs(val)<1e-30: zs.append(S)
    print("N",N,"total",tot,"proper vanishing subsums",zs)
# 3. golden block J_r = f_{r+1}+i f_r as 2-term sum of units of Q(i,sqrt5) with fixed coefficients
phi=(1+sp.sqrt(5))/2; phib=(1-sp.sqrt(5))/2; I=sp.I
a=(phi+I)/sp.sqrt(5); b=-(phib+I)/sp.sqrt(5)
ok=all(sp.simplify(a*phi**r+b*phib**r-(sp.fibonacci(r+1)+I*sp.fibonacci(r)))==0 for r in range(0,12))
print("J_r = a phi^r + b phibar^r for r<12:",ok)
# 4. gcd overlap of golden blocks with free exponents: |J_r|^2=f_{2r+1}; gcd(f_{2r+1},f_{2s+1})=f_gcd
from math import gcd
r_,s_=3*40+1,3*41+1  # 2r+1=243, 2s+1=249 -> gcd 3
print("gcd f_243,f_249 =",gcd(sp.fibonacci(243),sp.fibonacci(249)))
r_,s_=  (45-1)//2,(135-1)//2
print("2r+1=45,2s+1=135: log gcd/log f_45 =",float(sp.log(gcd(sp.fibonacci(45),sp.fibonacci(135)))/sp.log(sp.fibonacci(45))))
