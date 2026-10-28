import sympy as sp
m,Y=sp.symbols('m Y')
c = m**3-11*m**2+39*m-13
u = -(m+3)*c/(16*(m-3)**2*(m-1)**2); v = -(m+1)*c/(4*(m-3)**3*(m-1))
G = sp.Poly(sp.numer(sp.together((Y**2-Y+u)**2 + v*Y)), Y)
degs=[sp.degree(co,m) for co in G.all_coeffs()]
print("coefficient degrees in m:", degs, "max d =", max(degs))
D=sp.Poly(sp.discriminant(G),m); print("deg disc =", D.degree(), "; order at infinity = 6*d - deg =", 6*max(degs)-D.degree())
