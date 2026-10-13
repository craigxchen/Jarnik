# (1) k with c_1(F_G) = k sigma_1 on G(2,4) for the V_A foliation (exact, sympy):
#     along a general pencil of lines R_s (d_R = 1), deg F|_C = 2 d_R - deg prim(p_R ^ m'(0)).
# (2) d_R(beta) = Phi^* sigma_1 . beta from Plucker degrees of R(t) along Mobius test curves (mod p).
import sys, random, itertools
import sympy as sp
from sympy import Matrix, Rational as R, Poly, symbols
from star import point_of, quadric_through, LAB
from arr import PLANES, CHOSEN
s = symbols('s')
pts = [point_of([PLANES[i], PLANES[j], PLANES[k]]) for (i, j, k) in CHOSEN]
qb = quadric_through(pts)
def wedge(u, v):
    return [sp.expand(u[i]*v[j] - u[j]*v[i]) for i in range(len(u)) for j in range(i+1, len(u))]
def prim_deg(vec):
    vec = [sp.expand(c) for c in vec]
    nz = [c for c in vec if c != 0]
    g = nz[0]
    for c in nz[1:]: g = sp.gcd(g, c)
    return max(Poly(sp.cancel(c/g), s).degree() for c in nz)
def tangent_data(p, pp):
    # quadric through line span(p,pp) and the 6 points
    rows = [[sp.expand((u.T*Q*v)[0]) for Q in qb] for (u, v) in [(p, p), (p, pp), (pp, pp)]]
    Cm = Matrix(rows); cs = []
    for j in range(4):
        sub = Cm[:, [k for k in range(4) if k != j]]
        cs.append(sp.expand((-1)**j * sub.det(method='berkowitz')))
    Q = sum((cs[k]*qb[k] for k in range(4)), sp.zeros(4, 4))
    Q = Q.applyfunc(sp.expand)
    def other_line_pt(a, b):
        # w with R'_a = span(a, w): z in T_a Q, w = Q(z) b - 2 B(b,z) z
        tau = (a.T*Q)
        z = Matrix([tau[1], -tau[0], 0, 0])
        Qz = sp.expand((z.T*Q*z)[0]); Bbz = sp.expand((b.T*Q*z)[0])
        return (Qz*b - 2*Bbz*z).applyfunc(sp.expand)
    w = other_line_pt(p, pp); w2 = other_line_pt(pp, p)
    # q(t) = p + t w ; q2(t) = ((p+tw)^T Q w2) pp - t (w^T Q pp) w2
    pQw2 = sp.expand((p.T*Q*w2)[0]); wQw2 = sp.expand((w.T*Q*w2)[0]); wQpp = sp.expand((w.T*Q*pp)[0])
    q2_0 = pQw2*pp; q2_1 = wQw2*pp - wQpp*w2
    m0 = wedge(p, q2_0)
    m1 = [sp.expand(a + b) for a, b in zip(wedge(w, q2_0), wedge(p, q2_1))]
    return m0, m1
if __name__ == '__main__':
    random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
    rnd = lambda: R(random.randint(-20, 20), random.randint(1, 7))
    P0 = Matrix([rnd() for _ in range(4)]); a = Matrix([rnd() for _ in range(4)]); b = Matrix([rnd() for _ in range(4)])
    p = P0; pp = (a + s*b).applyfunc(sp.expand)
    m0, m1 = tangent_data(p, pp)
    dR = prim_deg(m0)
    W = wedge(m0, m1)
    dW = prim_deg(W)
    print('pencil: d_R =', dR, ' deg prim(p_R ^ tangent) =', dW, ' => k = deg F|_line =', 2*dR - dW)
