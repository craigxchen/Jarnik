"""Test of the even-M pinning theorem on the four-point Pell family
(integer_cotangent_lcm_height_target.md, Section 5): L=24,
X=(90UV+30V^2+3, 96UV, 160V^2+128UV+16), U+V sqrt5=(9+4 sqrt5)^n, n=1 mod 10.
We compute the primitive tuple, diameter/R^(1/3), and the cluster direction
theta (argument of the anchor point); the theorem predicts that 4*theta
converges (mod 2 pi) to the argument of a bounded Gaussian integer."""
import math
from itertools import combinations
from cot_lib import primitive_tuple, all_edge_lcm, gnorm, gsub, is_clique, gmul, gpow

def pell(n):
    U, V = 1, 0
    for _ in range(n):
        U, V = 9 * U + 20 * V, 4 * U + 9 * V
    return U, V

for n in (1, 11, 21, 31, 41):
    U, V = pell(n)
    L = 24
    X = sorted([90 * U * V + 30 * V * V + 3, 96 * U * V, 160 * V * V + 128 * U * V + 16])
    assert is_clique(X, L)
    rows = primitive_tuple(X, L)
    N = gnorm(rows[0])
    assert N == all_edge_lcm(X, L)
    R = math.sqrt(N)
    diam = max(math.sqrt(gnorm(gsub(a, b))) for a, b in combinations(rows, 2))
    th = math.atan2(rows[0][1], rows[0][0])
    # direction of z^4 for each point
    d4 = [(4 * math.atan2(z[1], z[0])) % (2 * math.pi) for z in rows]
    print('n=%d log10 N=%.1f diam/R^(1/3)=%.4f  4*theta (mod 2pi) per point=%s' %
          (n, math.log10(N), diam / R ** (1 / 3), ['%.12f' % t for t in d4]))
