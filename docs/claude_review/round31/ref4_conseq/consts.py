from fractions import Fraction as Fr
import itertools, math, random, sympy as sp
# (b) excess forced on a sample of k marks with E=0: delta >= 3(k-24)/(4k(k+3)) L
f=lambda k: Fr(3*(k-24),4*k*(k+3))
print("k=48:",f(48)," k=49:",f(49)," k=50:",f(50)," 1/f(49)=",float(1/f(49)))
best=max(range(25,200),key=f); print("argmax integer k:",best,f(best))
k=sp.symbols('k'); print("real critical point:",sp.solve(sp.diff(sp.Rational(3,4)*(k-24)/(k*(k+3)),k),k))
# Turan: complement K_k-free -> bad pairs >= C(m,2) - (1-1/(k-1)) m^2/2 ; fraction -> 1/(k-1)
for kk in (48,49):
    m=10**6; bad=Fr(math.comb(m,2))-Fr(kk-2,kk-1)*Fr(m*m,2)
    print(kk,"fraction at m=1e6:",float(bad/math.comb(m,2)),"limit 1/",kk-1)
# trade-off: product excess*fraction
print("excess*fraction best k:",max(range(25,200),key=lambda k:f(k)/(k-1)))
# (c) general sample size: Ebar >= [(3/4)(k-24) + k(m-k)/(4(m-1))]/(2k(k-1)) L
m_,k_=sp.symbols('m k',positive=True)
expr=(sp.Rational(3,4)*(k_-24)+k_*(m_-k_)/(4*(m_-1)))/(2*k_*(k_-1))
print("k=36 simplify minus target:",sp.simplify(expr.subs(k_,36)-(sp.Rational(1,140)-1/(8*(m_-1)))))
print("general minus (k-18)/(2k(k-1)) + 1/(8(m-1)):",sp.simplify(expr-((k_-18)/(2*k_*(k_-1))-1/(8*(m_-1)))))
g=lambda k: Fr(k-18,2*k*(k-1)); print("g(35),g(36),g(34),g(37):",g(35),g(36),g(34),g(37))
# hypergeometric second moment: exact enumeration check, min over x at x=1/2
def second_moment(m,k,s):
    tot=0; cnt=0
    for j in range(0,k+1):
        w=math.comb(s,j)*math.comb(m-s,k-j); tot+=w*Fr(j*2-k,2)**2; cnt+=w
    return tot/cnt
m,kk=50,36
vals={s:second_moment(m,kk,s) for s in range(0,m+1)}
smin=min(vals,key=vals.get); print("m=50,k=36 min second moment at s=",smin,vals[smin],"formula",Fr(kk*(m-kk),4*(m-1)))
for s in (0,10,25,37): 
    x=Fr(s,m); assert vals[s]==kk*x*(1-x)*Fr(m-kk,m-1)+(kk*x-Fr(kk,2))**2
print("variance formula exact for sampled s values")
# (c) pair separation -> sum (s_h-m/2)^2 <= mL/4 ; mean shared degree <= L/4 + L/(2 sqrt m)+O(L/m)
# random test of the identity sum_h s_h(m-s_h) = sum_{i<j} sep_ij
random.seed(3)
for trial in range(200):
    m=random.randint(4,12); L=random.randint(1,15)
    A=[[random.randint(0,1) for h in range(L)] for i in range(m)]
    lhs=sum(sum(A[i][h] for i in range(m))*(m-sum(A[i][h] for i in range(m))) for h in range(L))
    rhs=sum(sum(A[i][h]!=A[j][h] for h in range(L)) for i,j in itertools.combinations(range(m),2))
    assert lhs==rhs
print("separation identity ok")
# ratio missing/shared -> (1/140)/(1/4) = 1/35
print("limit missing fraction:",Fr(1,140)/Fr(1,4))
# finite-m version: lower bound on missing fraction with explicit error
def frac_lb(m):
    E=Fr(1,140)-Fr(1,8*(m-1))
    Gub=(Fr(m*m,4)+m*math.sqrt(m)/2+Fr(m,4))/(m*(m-1))   # per unit L, upper bound for mean shared degree (crude)
    return float(E)/Gub
for m in (36,100,1000,10**4,10**6): print("m=",m," missing fraction >=",frac_lb(m), " (1/35=",1/35,")")
