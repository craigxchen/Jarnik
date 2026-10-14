"""Independent certificate for the sdc7 (T=1, L=2, m=8) tangent octic: a real zero with large separation,
nonzero gradient and nonzero second fundamental form.  Exact rational arithmetic + explicit Lipschitz bounds."""
import json, mpmath as mp
from fractions import Fraction as Fr
import numpy as np
D = json.load(open('sdc7_forms.json'))
K = [k[1:] for k in D['0']['full']]          # relative coordinates k' (base 0), entries in {-1,0,1}
c = D['0']['c']; terms = [(ci, ki) for ci, ki in zip(c, K) if ci]
d = 7; m = 8; c1 = sum(abs(ci) for ci, _ in terms)
assert max(abs(x) for _, k in terms for x in k) <= 1
def lin(k, u): return sum(a*b for a, b in zip(k, u))
def Phi(u): return sum(ci * lin(k, u)**m for ci, k in terms)
def grad(u):
    g = [0]*d
    for ci, k in terms:
        l = lin(k, u); f = ci*m*l**(m-1)
        for i in range(d):
            if k[i]: g[i] += f*k[i]
    return g
def hess(u):
    H = [[0]*d for _ in range(d)]
    for ci, k in terms:
        l = lin(k, u); f = ci*m*(m-1)*l**(m-2)
        for i in range(d):
            if k[i]:
                for j in range(d):
                    if k[j]: H[i][j] += f*k[i]*k[j]
    return H
mp.mp.dps = 60
u0 = [0.50069288, -0.90210639, 0.75043512, 1.0, 0.24963353, -0.40126536, -0.65146347]
# Newton in mpmath along the gradient projected to {u_3 = 0}
u = [mp.mpf(x) for x in u0]
for it in range(30):
    g = grad(u); g[3] = 0
    gg = sum(x*x for x in g)
    t = -Phi(u)/gg
    u = [ui + t*gi for ui, gi in zip(u, g)]
print('Newton residual Phi/|c|_1 =', mp.nstr(Phi(u)/c1, 5))
ur = [Fr(int(mp.nint(x*10**30)), 10**30) for x in u]; ur[3] = Fr(1)
g = grad([float(x) for x in ur]); g[3] = 0.0
gr = [Fr(int(round(x/max(abs(y) for y in g)*10**6)), 10**6) for x in g]
h = Fr(1, 10**22)
pm = Phi([a - h*b for a, b in zip(ur, gr)]); pp = Phi([a + h*b for a, b in zip(ur, gr)])
print('sign change of Phi on the segment u_r +- h g (h = 1e-22):', (pm < 0) != (pp < 0), float(pm), float(pp))
# |u*|_inf = 1 exactly (u_3 = 1, g_3 = 0), other coordinates
print('max |u_x| over x != 3:', max(abs(float(x)) for i, x in enumerate(ur) if i != 3))
z = [Fr(0)] + ur
gaps = sorted(abs(z[i]-z[j]) for i in range(8) for j in range(i+1, 8))
gn = max(abs(x) for x in gr)
print('sep(u_r) = %.9f ; sep(u*) >= %.9f' % (float(gaps[0]), float(gaps[0] - 2*h*gn)))
# gradient at u*: |grad(u*) - grad(u_r)|_2 <= |u*-u_r|_2 * sup ||Hess||_F on the segment
# ||Hess||_F <= sum |c_k| m(m-1) |<k,x>|^{m-2} |k|_2^2 <= c1*56*(7*1.0000001)^6*7
dist = float(h) * float(sum(x*x for x in gr))**0.5
HB = c1*56*(7*1.0000001)**6*7
G = grad(ur); Gn = float(sum(x*x for x in G))**0.5
print('|grad Phi(u_r)|_2 = %.4e ; perturbation <= %.3e  => grad Phi(u*) != 0: %s' % (Gn, HB*dist, Gn > HB*dist))
# restricted Hessian on grad^perp at u_r, exact; third-derivative bound for the perturbation
H = np.array([[float(x) for x in row] for row in hess(ur)]); Gv = np.array([float(x) for x in G])
nG = Gv/np.linalg.norm(Gv); P = np.eye(d) - np.outer(nG, nG)
PHP = P @ H @ P
TB = c1*56*6*(7*1.0000001)**5*7**1.5     # bound on ||D^3 Phi||: sum |c| m(m-1)(m-2)|<k,x>|^{m-3}|k|_2^3
# variation of P H P: ||P*H*P* - P H P|| <= ||H*-H|| + 2||H|| ||P*-P||, ||P*-P|| <= 2|G*-G|/|G|
Hn = np.linalg.norm(H); var = TB*dist + 2*Hn*(2*HB*dist/Gn)
print('||P H P||_F at u_r = %.4e ; variation bound <= %.3e  => second fundamental form nonzero: %s' % (np.linalg.norm(PHP), var, np.linalg.norm(PHP) > var))
print('positions (0,u*) sorted:', sorted([0.0] + [float(x) for x in ur]))
