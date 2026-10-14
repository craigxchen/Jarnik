"""Joint system of the 8 base-point octics of sdc7 (T=1, L=2, m=8)."""
import json, random, flint, numpy as np
D = json.load(open('sdc7_forms.json'))
n = 8
F = [(np.array(D[str(b)]['full'], dtype=np.int64), np.array(D[str(b)]['c'], dtype=object)) for b in range(n)]
def phi_exact(b, z):   # z integer vector (python ints), exact
    K, c = F[b]
    lin = K.astype(object).dot(np.array(z, dtype=object))
    return int(sum(ci * li**8 for ci, li in zip(c, lin)))
random.seed(1)
pts = [[0] + [random.randint(-50, 50) for _ in range(n-1)] for _ in range(40)]
Mx = flint.fmpz_mat([[phi_exact(b, z) for b in range(n)] for z in pts])
print('exact rank of the 8 base-point octics (evaluated at 40 integer points):', Mx.rank())
# the author's sdc7 zero (rounded), relative to base 0: z = (0, u*)
u = [0.50069, -0.90211, 0.75044, 1, 0.24963, -0.40127, -0.65146]
z = np.array([0.0] + u)
for b in range(n):
    K, c = F[b]; lin = K @ z
    val = float(np.dot(c.astype(float), lin**8)) / (1008 * 7.0**8)
    print('  base %d: Phi_b(z*)/(|c|_1 d^8) = %+.3e' % (b, val))
