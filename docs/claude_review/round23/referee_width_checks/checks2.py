import sympy as sp, itertools
from fractions import Fraction
x=sp.symbols('x1:6')
def supp(expr, shift):
    P=sp.Poly(sp.expand(expr),*x)
    return [(tuple(e_-s for e_,s in zip(e,shift)),c) for e,c in P.terms()]
def order_at_e(terms, maxd=15):
    # ord >= t iff sum c_s P(s)=0 for all monomials P of deg < t (exact integers)
    k=len(terms[0][0])
    for d in range(maxd+1):
        for e in itertools.product(range(d+1),repeat=k):
            if sum(e)!=d: continue
            val=sum(int(c)*sp.prod([s[j]**e[j] for j in range(k)]) for s,c in terms)
            if val!=0: return d
    return None
V=sp.prod([x[i]-x[j] for i in range(5) for j in range(i+1,5)])
tV=supp(V,(2,)*5)
print("V: max sum_+ over support =",max(sum(max(a,0) for a in s) for s,_ in tV),"; all sum=0:",all(sum(s)==0 for s,_ in tV))
print("V: ord_e =",order_at_e(tV,11))
V3=(x[0]-x[1])*(x[0]-x[2])*(x[1]-x[2])
t3=supp(V3,(1,1,1,0,0))
print("V123: max sum_+ =",max(sum(max(a,0) for a in s) for s,_ in t3),"; ord_e =",order_at_e(t3,4))
# widths of V support over cuts
for r in (1,2):
    ws=set()
    for A in itertools.combinations(range(5),r):
        vals=[sum(s[j] for j in A) for s,_ in tV]; ws.add(max(vals)-min(vals))
    print("V widths |A|=%d:"%r,ws)
# template computation: peel of slice bounds for S_m sliced by x_1
def peel(vals):
    v=sorted(vals,reverse=True); return max(v[i]+i for i in range(len(v)))
print("m, peel with Cbal(4) floor(3/2(2m-|c|)), peel with ThmB(4) floor(7/4(2m-|c|)), 7m/2, 15m/4, floor(10m/3)")
for m in (1,2,3,4,6,8,12,24,48):
    cs=range(-m,m+1)
    p1=peel([ (3*(2*m-abs(c)))//2 for c in cs]); p2=peel([(7*(2*m-abs(c)))//4 for c in cs])
    print(m,p1,p2,Fraction(7*m,2),Fraction(15*m,4),10*m//3)
