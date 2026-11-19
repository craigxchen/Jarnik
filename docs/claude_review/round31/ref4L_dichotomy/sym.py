# Symbolic checks: interpolation identity, binary Vandermonde, Gram identity, profile formulas.
from sympy import symbols, Matrix, expand, simplify, factor, prod, Rational
from itertools import combinations
xs = symbols('x0:5'); ys = symbols('y0:5')
V = [(xs[i], ys[i]) for i in range(5)]
def det(u, v): return u[0]*v[1]-u[1]*v[0]
a, b, c = symbols('a b c')
Q = lambda v: a*v[0]**2 + b*v[0]*v[1] + c*v[1]**2
# interpolation identity sum_i Q(V_i) V_i / prod_{j!=i} m_ij, m_ij = [V_i,V_j]
for coord in range(2):
    s = sum(Q(V[i])*V[i][coord] / prod([det(V[i], V[j]) for j in range(5) if j != i]) for i in range(5))
    assert simplify(s) == 0
print("interpolation identity: OK (generic symbolic, both coordinates)")
# binary Vandermonde
M = Matrix([[V[i][0]**2, V[i][0]*V[i][1], V[i][1]**2] for i in range(3)])
print("Vandermonde det / (m01 m02 m12) =", factor(M.det()/(det(V[0],V[1])*det(V[0],V[2])*det(V[1],V[2]))))
# Gram: for symmetric E, (E(Vi)E(Vj) - E(Vi,Vj)^2) = m_ij^2 det E
e1, e2, e3 = symbols('e1 e2 e3')
Em = Matrix([[e1, e2/2], [e2/2, e3]])
vi = Matrix(V[0]); vj = Matrix(V[1])
lhs = (vi.T*Em*vi)[0]*(vj.T*Em*vj)[0] - ((vi.T*Em*vj)[0])**2
assert expand(lhs - det(V[0], V[1])**2*Em.det()) == 0
print("Gram identity: OK")
# V_k = a V_i + b V_j => 2E(Vi,Vj) = (E(Vk) - a^2 E(Vi) - b^2 E(Vj))/(ab)
vk = Matrix(V[2])
aa = det(V[2], V[1])/det(V[0], V[1]); bb = det(V[0], V[2])/det(V[0], V[1])
assert simplify(aa*vi + bb*vj - vk) == Matrix([0, 0])
E = lambda u, w: (u.T*Em*w)[0]
assert simplify(2*E(vi, vj) - (E(vk, vk) - aa**2*E(vi, vi) - bb**2*E(vj, vj))/(aa*bb)) == 0
print("polarization formula: OK")
# profile formulas n_i = D B_i/A_i^2, g_ij = D J_ij/(A_iA_j), prod_{j!=i} g_ij = D^2 B_i/A_i^3, both orientations
edges = list(combinations(range(5), 2)); J = {e: symbols('J%d%d' % e) for e in edges}
for F in [set(), {(0,1),(1,2),(2,3),(3,4),(0,4)}]:
    S = [{e for e in edges if (i in e) != (e in F)} for i in range(5)]
    D = prod([J[e] for e in F]); A = [prod([J[e] for e in F if i in e]) for i in range(5)]
    B = [prod([J[e] for e in edges if i in e]) for i in range(5)]
    n = [prod([J[e] for e in S[i]]) for i in range(5)]
    for i in range(5):
        assert simplify(n[i] - D*B[i]/A[i]**2) == 0
        g = [prod([J[e] for e in S[i] & S[j]]) for j in range(5)]
        for j in range(5):
            if j != i:
                e = tuple(sorted((i, j)))
                assert simplify(g[j] - D*J[e]/(A[i]*A[j])) == 0
        assert simplify(prod([g[j] for j in range(5) if j != i]) - D**2*B[i]/A[i]**3) == 0
    marks = {e: sum(1 for i in range(5) if e in S[i]) for e in edges}
    print("orientation |F|=%d: profile formulas OK; |S_i| =" % len(F), [len(s) for s in S],
          "|S_i cap S_j| =", sorted({len(S[i]&S[j]) for i,j in edges}), "marked per edge:", sorted(set(marks.values())))
