"""T16: base-coordinate rigidity (R) on RANDOM balanced sign-column profiles (no symmetry), by two IPs:
  h = max sigma over A;   feasibility of {k in A, sigma(k) = h, k + e_x - e_y in A for SOME ordered pair x != y}
(pair chosen by binary variables z_x (+1 at x), u_y (-1 at y)).  (R) <=> infeasible.
Profile: r patterns, each a uniformly random subset of the d = n-1 non-base points of size floor(d/2) or ceil(d/2)
(balanced at the base up to O(1/sqrt r) fluctuations).  Usage: t16_random_rigidity.py n r L seed timelimit"""
import sys, time, random
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
n, r, L, seed, tl = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), float(sys.argv[5])
d = n - 1
rng = random.Random(seed)
Pm = np.zeros((r, d))
for p in range(r):
    for x in rng.sample(range(d), d // 2 + rng.randint(0, d % 2)): Pm[p, x] = 1
t0 = time.time()
res = milp(-np.ones(d), constraints=[LinearConstraint(Pm, -L, L)], integrality=np.ones(d), bounds=Bounds(-1, 1),
           options={'time_limit': tl})
h = round(-res.fun)
print("random profile n=%d r=%d L=%d seed=%d: h = %d (status %d, bound %s) (%.0fs)" % (n, r, L, seed, h, res.status,
      getattr(res, 'mip_dual_bound', None), time.time() - t0), flush=True)
if res.status != 0: sys.exit(0)
# variables: k (d, in [-1,1]), z (d, binary), u (d, binary)
N = 3 * d
I = np.eye(d); Z = np.zeros((d, d)); Zr = np.zeros((r, d))
cons = [LinearConstraint(np.hstack([Pm, Zr, Zr]), -L, L),                       # k in A
        LinearConstraint(np.hstack([Pm, Pm, -Pm]), -L, L),                      # k + z - u in A
        LinearConstraint(np.hstack([np.ones((1, d)), np.zeros((1, d)), np.zeros((1, d))]), h, h),   # sigma(k) = h
        LinearConstraint(np.hstack([np.zeros((1, d)), np.ones((1, d)), np.zeros((1, d))]), 1, 1),   # one +1
        LinearConstraint(np.hstack([np.zeros((1, d)), np.zeros((1, d)), np.ones((1, d))]), 1, 1),   # one -1
        LinearConstraint(np.hstack([I, I, -I]), -1, 1),                         # k + z - u in [-1,1]
        LinearConstraint(np.hstack([Z, I, I]), 0, 1)]                           # x != y
lb = np.concatenate([-np.ones(d), np.zeros(2 * d)]); ub = np.ones(N)
res2 = milp(np.zeros(N), constraints=cons, integrality=np.ones(N), bounds=Bounds(lb, ub), options={'time_limit': tl})
print("   top-level pair move exists: %s (status %d) -> (R) %s (%.0fs)" % (res2.status == 0, res2.status,
      "FAILS" if res2.status == 0 else ("HOLDS" if res2.status == 2 else "undecided"), time.time() - t0), flush=True)
