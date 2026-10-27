# Judge check J1 (exact, sympy rationals): toy refutation of "all probabilistic forms equivalent".
# Psi_S(u) := min_{P deg<m} sum_{k in S} (l_u(k)^m - P(k))^2   (sum form; mean form = /|S|)
import itertools, sympy as sp
u1,u2=sp.symbols('u1 u2')
def monos(d,m):
    return [a for a in itertools.product(range(m),repeat=d) if sum(a)<m]
def psi(S,m,u):
    S=list(S); d=len(S[0])
    M=sp.Matrix([[sp.prod([sp.Integer(k[i])**a[i] for i in range(d)]) for a in monos(d,m)] for k in S])
    g=sp.Matrix([sum(sp.Integer(k[i])*u[i] for i in range(d))**m for k in S])
    # projection residual via column space basis
    C=M.columnspace(); B=sp.Matrix.hstack(*C)
    coef=(B.T*B).LUsolve(B.T*g)
    r=g-B*coef
    return sp.factor(sp.simplify((r.T*r)[0]))
cube=list(itertools.product([-1,0,1],repeat=2))
A=[k for k in cube if abs(k[0]+k[1])<=1]
print("|A|=",len(A))
PA=psi(A,3,(u1,u2)); print("Psi_A(u) symbolic (sum form):",PA)
print("Psi_A(1/2,1) sum form =",PA.subs({u1:sp.Rational(1,2),u2:1}),"; mean form =",PA.subs({u1:sp.Rational(1,2),u2:1})/len(A))
B=[k for k in cube if k not in [(1,1),(-1,0)]]
print("|B|=",len(B)," Psi_B(1/2,1)=",psi(B,3,(sp.Rational(1,2),1)))
A6=[k for k in A if k!=(-1,0)]
print("|A6|=",len(A6)," Psi_A6(1/2,1)=",psi(A6,3,(sp.Rational(1,2),1)))
print("cube m=4 Psi(u) =",psi(cube,4,(u1,u2)))
cm=[k for k in cube if k!=(1,1)]
print("cube minus point m=4 Psi(u) =",psi(cm,4,(u1,u2)))
