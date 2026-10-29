# Family generalising the C=7/2 curve: conic u(u+a+b) = k w(w+c+d), k = ab/(cd);
# five divisors D1=(u+a)(u+b), D2=(u+a)(u+a+b), D3=(u+b)(u+a+b), D4=k(w+c)(w+c+d), D5=k(w+d)(w+c+d).
# Exact symbolic Jacobian rank of (a,b,c,d,t) -> (D_j/D_1)_{j=2..5} at rational points.
import sympy as sp, random
a,b,c,d,t=sp.symbols('a b c d t')
k=a*b/(c*d)
u=(k*t*(c+d)-a-b)/(1-k*t**2); w=t*u
assert sp.simplify(u*(u+a+b)-k*w*(w+c+d))==0
D=[(u+a)*(u+b),(u+a)*(u+a+b),(u+b)*(u+a+b),k*(w+c)*(w+c+d),k*(w+d)*(w+c+d)]
# identity check of the second factorisation (u+a)(u+b) = k(w+c)(w+d)
print("second factorisation identity:", sp.simplify((u+a)*(u+b)-k*(w+c)*(w+d))==0)
F=[sp.simplify(D[j]/D[0]) for j in range(1,5)]
vars_=[a,b,c,d,t]
J=sp.Matrix([[sp.diff(f,v) for v in vars_] for f in F])
random.seed(1)
for trial in range(4):
    pt={a:sp.Rational(random.randint(1,9),random.randint(1,5)),b:sp.Rational(random.randint(1,9),random.randint(1,5)),
        c:sp.Rational(random.randint(1,9),random.randint(1,5)),d:sp.Rational(random.randint(1,9),random.randint(1,5)),
        t:sp.Rational(random.randint(1,20),random.randint(1,20))}
    Jn=J.subs(pt)
    print(trial, {str(k_):str(v) for k_,v in pt.items()}, "rank", Jn.rank())
