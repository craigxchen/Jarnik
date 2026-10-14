"""sdc7 (T=1, L=2, m=8): the 8 base-point systems have the same octic Phi (exact), but different exponential sums.
Second-order terms Psi_b(u) = sum_k c^b_k <k,u>^10; joint zeros of all F_{c_b} near e need Phi = O(t^2) and
Psi_b - Psi_0 = O(t^2).  Exact: is Psi_b - Psi_0 in the ideal (Phi) in degree 10? Evaluate at the author's zero."""
import json, random, flint, numpy as np
from fractions import Fraction
D = json.load(open('sdc7_forms.json')); n = 8
F = []
for b in range(n):
    K = np.array(D[str(b)]['full'], dtype=object); c = np.array(D[str(b)]['c'], dtype=object); F.append((K, c))
def mom(b, z, p):
    K, c = F[b]; lin = K.dot(np.array(z, dtype=object))
    return sum(ci * li**p for ci, li in zip(c, lin))
random.seed(3)
z0 = [0] + [random.randint(-20, 20) for _ in range(n-1)]
sgn = [1 if mom(b, z0, 8) == mom(0, z0, 8) else (-1 if mom(b, z0, 8) == -mom(0, z0, 8) else None) for b in range(n)]
print('normalisation signs (Phi_b = sgn*Phi_0 exactly at a random point):', sgn)
# exact identity check at 30 random points
ok = all(sgn[b]*mom(b, z, 8) == mom(0, z, 8) for b in range(n) for z in ([0]+[random.randint(-30,30) for _ in range(n-1)] for _ in range(5)))
print('Phi_b == sgn_b Phi_0 at 40 random integer points:', ok)
# rank of {Psi_b - Psi_0} modulo the ideal (Phi) in degree 10: compare rank of [Psi_b - Psi_0 ; Phi*monomials of degree 2]
pts = [[0] + [random.randint(-40, 40) for _ in range(n-1)] for _ in range(80)]
quad = [(i, j) for i in range(1, n) for j in range(i, n)]   # z_0 = 0 coordinates
rowsI = []; rowsD = []
for z in pts:
    ph = mom(0, z, 8)
    rowsI.append([ph * z[i] * z[j] for (i, j) in quad])
    rowsD.append([sgn[b]*mom(b, z, 10) - mom(0, z, 10) for b in range(1, n)])
MI = flint.fmpz_mat(rowsI); MID = flint.fmpz_mat([a + bb for a, bb in zip(rowsI, rowsD)]); MD = flint.fmpz_mat(rowsD)
print('rank Phi*quadrics:', MI.rank(), ' rank {Psi_b - Psi_0}:', MD.rank(), ' rank together:', MID.rank())
print(' => dim of span{Psi_b - Psi_0} modulo (Phi)_10 =', MID.rank() - MI.rank())
# evaluate at the author's (rounded) zero
u = [0.50069288, -0.90210639, 0.75043512, 1.0, 0.24963353, -0.40126536, -0.65146347]
z = np.array([0.0] + u)
def momf(b, z, p):
    K, c = F[b]; lin = K.astype(float) @ z
    return float(np.dot(c.astype(float), lin**p))
norm8 = 1008 * 7.0**8; norm10 = 1008 * 7.0**10
print('Phi(u*)/(|c|_1 d^8) = %.2e' % (momf(0, z, 8)/norm8))
for b in range(1, n):
    print('  b=%d: (Psi_b - Psi_0)(u*)/(|c|_1 d^10) = %+.4e   [Psi_0(u*)/(|c|_1 d^10) = %+.4e]' % (b, (sgn[b]*momf(b, z, 10) - momf(0, z, 10))/norm10, momf(0, z, 10)/norm10))
