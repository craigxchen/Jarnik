# exact symbolic checks of the reformulation identities (sympy)
import sympy as sp, itertools
s=sp.symbols('s',real=True)
# (1) Moebius converse identity: z_a = -i det(zhat_a, g); Im(z_a conj z_b) = Im(g1 conj g2) * det(zhat_a, zhat_b)
p=sp.symbols('p0:4',real=True); q=sp.symbols('q0:4',real=True)
zha=(p[0]*s+p[1], p[2]*s+p[3]); zhb=(q[0]*s+q[1], q[2]*s+q[3])
G=[sp.symbols(f'g{k}r0:5',real=True) for k in (1,2)]; H=[sp.symbols(f'g{k}i0:5',real=True) for k in (1,2)]
g=[sum((G[k][d]+sp.I*H[k][d])*s**d for d in range(5)) for k in range(2)]
gb=[sum((G[k][d]-sp.I*H[k][d])*s**d for d in range(5)) for k in range(2)]
za=-sp.I*(zha[0]*g[1]-zha[1]*g[0]); zb=-sp.I*(zhb[0]*g[1]-zhb[1]*g[0])
zbb=sp.I*(zhb[0]*gb[1]-zhb[1]*gb[0])   # conj(z_b) coefficientwise (s real)
lhs=sp.expand(za*zbb); im_lhs=sp.expand((lhs-sp.expand(sp.conjugate(lhs)))/(2*sp.I))
im_g=sp.expand((g[0]*gb[1]-gb[0]*g[1])/(2*sp.I))
det_ab=sp.expand(zha[0]*zhb[1]-zha[1]*zhb[0])
print('(1) Moebius identity holds:', sp.simplify(im_lhs-sp.expand(im_g*det_ab))==0)
# (2) half-angle identity C_a = z_a * Omegabar / zbar_a with random numeric places (exact Gaussian rationals)
import random
random.seed(1)
FULL=frozenset(range(1,6))
small=[frozenset([a]) for a in range(1,6)]+[frozenset(c) for c in itertools.combinations(range(1,6),2)]
rho={A: sp.Rational(random.randint(-9,9),random.randint(1,5))+sp.I*sp.Rational(random.choice([-1,1])*random.randint(1,9),random.randint(1,5)) for A in small}
def place(A):
    return rho[A] if A in rho else sp.conjugate(rho[FULL-A])
allA=[frozenset(c) for k in range(1,5) for c in itertools.combinations(range(1,6),k)]
ok=True
Om=sp.prod([s-rho[A] for A in small]); Omb=sp.prod([s-sp.conjugate(rho[A]) for A in small])
for a in range(1,6):
    Ca=sp.prod([s-place(A) for A in allA if a in A])
    za=(s-rho[frozenset([a])])*sp.prod([s-rho[frozenset([a,x])] for x in range(1,6) if x!=a])
    zab=sp.expand(sp.conjugate(za).subs(sp.conjugate(s),s))
    ok&= sp.expand(Ca*zab - za*Omb)==0
print('(2) C_a * zbar_a == z_a * Omegabar for all a:', ok)
# (3) column norm: C_a * Cbar_a = n(s) independent of a
ns=[sp.expand(sp.prod([s-place(A) for A in allA if a in A])*sp.prod([s-sp.conjugate(place(A)) for A in allA if a in A])) for a in range(1,6)]
print('(3) C_a Cbar_a independent of a:', all(sp.expand(n-ns[0])==0 for n in ns))
