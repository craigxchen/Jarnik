"""Classify pentagon-signed first and second complex star moments."""
from sympy import Matrix, Symbol, symbols

# u = positive cycle, with u4 set to zero by a common translation.
u=list(symbols('a:4'))+[0]
d=Symbol('d')
V=u[:4]+[d]
v=[u[j]+u[(j+1)%5]-u[(j-2)%5]+d for j in range(5)]
q=[(v[j]**2-u[j]**2-u[(j+1)%5]**2+u[(j-2)%5]**2)/2 for j in range(5)]
F=Matrix([q[j]-q[0] for j in range(1,5)])
J=F.jacobian(V)
for b,c,t in [(0,0,1),(1,0,0),(-1,2,2)]:
    z=Matrix([0,b,c,b,t])
    fix=dict(zip(V,z))
    assert F.subs(fix)==Matrix.zeros(4,1)
    A=J.subs(fix)
    nonzero_minor=A.extract(A.T.rref()[1],A.rref()[1]).det()
    print('representative',list(z),'Jacobian minor',nonzero_minor)
    for vec in A.nullspace():
        # Vertex reflection i->-i flips u_i to u_(-i-1).
        # Restore gauge u4=0 by subtracting u0 from every edge coordinate.
        a0,a1,a2,a3,dv=vec
        reflected=Matrix([-a0,a3-a0,a2-a0,a1-a0,dv])
        assert reflected==vec
    if A.rank()==3:
        k=A.nullspace()[1]
        obstruction=(A.T.nullspace()[0].T*F.subs(dict(zip(V,k))))[0]
        print('quadratic obstruction at double point',obstruction)
        assert obstruction==2
