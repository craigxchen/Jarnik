# Independent symbolic checks of IV.2 (generic x with s4=0 imposed via elimination where needed)
from sympy import *
from itertools import combinations
P,q=symbols('P q')
X=symbols('x1:6')
s=[sum(prod(c) for c in combinations(X,k)) for k in range(6)]
s1,s2,s3,s4,s5=s[1:]
b=lambda i,j: X[i]*X[j]*(s1-X[i]-X[j])
h=[X[i]**3-s1*X[i]**2+s2*X[i] for i in range(5)]
# identity 1
for i,j in combinations(range(5),2):
    assert expand(X[i]*h[j]-X[j]*h[i]-(X[i]-X[j])*b(i,j))==0
print("x_i h_j - x_j h_i = (x_i-x_j) b_ij : generic OK")
# shell identities generic (before s4=0)
i=0
others=[j for j in range(5) if j!=i]
B=[b(i,j) for j in others]
e1=sum(B); e2=sum(B[a]*B[c] for a,c in combinations(range(4),2))
e3=sum(B[a]*B[c]*B[d] for a,c,d in combinations(range(4),3)); e4=prod(B)
x=X[0]
assert expand(e1-2*h[0])==0
d2=expand(h[0]**2-e2-(s4*x**2-(s1*s4+3*s5)*x+s1*s5))
print("generic second identity residual:",d2)
# Now norms with symbolic l,b0 such that l^2=-s1 s5, b0=3 s5/2
l2=-s1*s5; b0=Rational(3,2)*s5
# |V|^2 = P^4 + 2h P^3 q + (h^2+2 b0 x + l^2) P^2 q^2 + (2 b0 x h + 2 l^2 h) P q^3 + (b0^2 x^2 + l^2 h^2) q^4
c2=h[0]**2+2*b0*x+l2
print("P^2q^2 coefficient minus e2, generic:",factor(expand(c2-e2)))
c3=2*b0*x*h[0]+2*l2*h[0]
r3=expand(c3-e3)
claim=s5*x*(s2*x-s1*s2-s3)
diff=expand(r3-claim)
# reduce diff modulo s4=0: substitute x5 solving s4=0 is rational: s4 linear in x5
x5=X[4]
sol=solve(Eq(s4,0),x5)[0]
print("P q^3 coefficient minus claim, mod s4=0:",simplify(together(diff.subs(x5,sol))))
print("diff mod s4 (as poly in s4?) :",factor(diff))
