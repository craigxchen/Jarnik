# Independent check of the five-cycle relabelling case (III.7): u_{i,i+j}=A_j z^i, v_{i,i+j}=B_j z^i.
from sympy import symbols, groebner, expand, Rational, simplify
z,A1,B1,A2,B2,t=symbols('z A1 B1 A2 B2 t')
Phi=z**4+z**3+z**2+z+1
# star sums of u^k-v^k at mark i: bracket must vanish for k=1..4 (k=5 automatically i-independent).
eqs=[]
for k in range(1,5):
    # star at mark 0: edges {0,1}(A1 z^0),{4,0}(A1 z^4),{0,2}(A2 z^0),{3,0}(A2 z^3)
    br=(A1**k-B1**k)*(1+z**(4*k%5))+(A2**k-B2**k)*(1+z**(3*k%5))
    eqs.append(expand(br))
# strictness: A1!=B1 (delta1!=0) -> Rabinowitsch with t; also include Phi
G=groebner(eqs+[Phi, t*(A1-B1)-1, A1-B1-2],A1,B1,A2,B2,t,z,order='lex')
print("normalized (A1-B1=2) Groebner basis:",G.exprs[:3], "... len",len(G.exprs))
