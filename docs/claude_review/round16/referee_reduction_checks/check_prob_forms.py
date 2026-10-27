# Toy counterexamples to "all probabilistic versions are the same statement up to e^{O(d)}" (reduction.md Cor 2.2 / summary 1).
# d = 2 (n = 3, base point u_3 = 0), m = 3, slab box A = {k' in {-1,0,1}^2 : |k_1 + k_2| <= 1} (one sign column, S = {1,2}, L = 1;
# equivalently |<k,beta>| <= 2 with beta = (1,1,-1), k_3 = -k_1-k_2).
import itertools, sys
from fractions import Fraction
sys.path.insert(0, '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round16/referee_reduction_checks')
from rlib import psi_sigma

cube = list(itertools.product((-1, 0, 1), repeat=2))
A = [k for k in cube if abs(k[0] + k[1]) <= 1]
B = [k for k in cube if k not in [(1, 1), (-1, 0)]]      # another 7-point subset of the cube
m = 3
u = [Fraction(1, 2), Fraction(1)]                         # |u|_inf = 1, pairwise gaps 1/2, 1/2, 1 (no collision)
print('|A| =', len(A), ' |B| =', len(B))
print('Psi^Sigma_A(u) =', psi_sigma(A, m, u), ' (> 0: Cheb_A(u) > 0, MIN_inf holds at u for small eps_0)')
print('Psi^Sigma_B(u) =', psi_sigma(B, m, u), ' (= 0: some P with l_u^3 = P on B, cube mass 7/9 = mu(A))')
# TC(A,3) on a grid: is it inside the braid arrangement u1*u2*(u1-u2) = 0 ?
bad = []
for a in range(-6, 7):
    for b in range(-6, 7):
        if (a, b) == (0, 0): continue
        if psi_sigma(A, m, [a, b]) == 0 and a * b * (a - b) != 0:
            bad.append((a, b))
print('grid points of TC(A,3) off the braid arrangement:', bad)
# AC_p failure: a subset of A of size 6 (mass 6/7 under mu_A) on which u is a tangent-cone point
for drop in A:
    Bp = [k for k in A if k != drop]
    if psi_sigma(Bp, m, u) == 0:
        print('subset A minus', drop, ': Psi = 0, so SB_{mu_A}(u;0) >= 6/7 (AC_p fails for every p > 1/7, every delta)')
        break
# general d: cube minus one point at m = rho (top degree): all u become TC points
C2 = cube
C2m = [k for k in cube if k != (0, 0)]
print('cube d=2, m=4: Psi_cube(u) =', psi_sigma(C2, 4, u), '; cube minus (0,0): Psi =', psi_sigma(C2m, 4, u))
