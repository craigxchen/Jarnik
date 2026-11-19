# Symbolic sanity check of the algebra in the rigidity proof (Gauss-lemma factorisation and the
# rank-one exclusion): for random primitive integer V_a and integral L, DE=[V_a,X]L(X) has
# content(DE)=content(L); if additionally L(V_j)=0 for primitive V_j then L=m[V_j,X].
import random, math
from sympy import symbols, Poly, gcd as sgcd
X,Y=symbols('X Y')
def prim():
    while True:
        a,b=random.randint(-10**6,10**6),random.randint(-10**6,10**6)
        if math.gcd(a,b)==1: return a,b
bad=0
for trial in range(2000):
    a=prim(); l=(random.randint(-50,50),random.randint(-50,50))
    if l==(0,0): continue
    Va=a[1]*X-a[0]*Y            # [V_a,(X,Y)] = a0*Y - a1*X up to sign
    Lf=l[0]*X+l[1]*Y
    DE=Poly((Va*Lf).expand(),X,Y)
    cont=math.gcd(*[int(c) for c in DE.coeffs()])
    if cont!=math.gcd(*l): bad+=1
print('Gauss-lemma content test failures:',bad)
