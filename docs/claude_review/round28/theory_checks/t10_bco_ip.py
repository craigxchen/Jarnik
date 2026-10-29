"""T10: base-coordinate obstruction at larger Sylvester sizes by integer programming (scipy.milp / HiGHS).
Sylvester n = 2^k, T = 1 box, slabs |k(P_u)| <= L (P_u = side omitting the base).  Computes
  h      = max sigma over A                                   (sigma = sum of the d = n-1 cube coordinates = -k_n)
  h_pair = max sigma over {k : k in A, k + e_x - e_y in A}   (one pair suffices: GL(k,2) fixing the base is
                                                              2-transitive on the other points)
  z2     = max |S1|+|S2|, S1,S2 disjoint 0/1 sets with |S_i cap P_u| <= L   (box-zonotope certificate of (N))
(R) <=> h_pair < h.  If (R) and z2 >= 2h then reach(A, L_pair) <= 2h <= zeta_pair(A) (Prop BCO).
Usage: t10_bco_ip.py k L [timelimit]"""
import sys, itertools, time
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
k, L = int(sys.argv[1]), int(sys.argv[2]); tl = float(sys.argv[3]) if len(sys.argv) > 3 else 600
n = 2 ** k; d = n - 1; base = n - 1
pts = list(range(n))
def dot(u, x): return bin(u & x).count('1') % 2
Pm = []
for u in range(1, n):
    S = [x for x in pts if dot(u, x) == 1]
    P = S if base not in S else [x for x in pts if x not in S]
    Pm.append([1 if x in P else 0 for x in range(d)])
Pm = np.array(Pm, dtype=float)       # (n-1) x d
c = -np.ones(d)
def solve_h():
    res = milp(c, constraints=[LinearConstraint(Pm, -L, L)], integrality=np.ones(d), bounds=Bounds(-1, 1),
               options={'time_limit': tl, 'disp': False})
    return res
def solve_pair(x=0, y=1):
    sh = Pm[:, x] - Pm[:, y]
    lb = np.full(d, -1.0); ub = np.full(d, 1.0); ub[x] = 0; lb[y] = 0
    res = milp(c, constraints=[LinearConstraint(Pm, -L, L), LinearConstraint(Pm, -L - sh, L - sh)],
               integrality=np.ones(d), bounds=Bounds(lb, ub), options={'time_limit': tl, 'disp': False})
    return res
def solve_z2():
    # variables s1 (d), s2 (d) in {0,1}; s1+s2 <= 1; Pm s_i <= L
    A1 = np.hstack([Pm, np.zeros_like(Pm)]); A2 = np.hstack([np.zeros_like(Pm), Pm])
    Dj = np.hstack([np.eye(d), np.eye(d)])
    res = milp(-np.ones(2 * d), constraints=[LinearConstraint(A1, -np.inf, L), LinearConstraint(A2, -np.inf, L),
               LinearConstraint(Dj, -np.inf, 1)], integrality=np.ones(2 * d), bounds=Bounds(0, 1),
               options={'time_limit': tl, 'disp': False})
    return res
t0 = time.time()
r = solve_h(); h = -r.fun; hb = getattr(r, 'mip_dual_bound', None)
print("Sylvester n=%d d=%d L=%d |P|=%d  LP-type bound 2L(1-1/n)=%.2f" % (n, d, L, n // 2, 2 * L * (1 - 1 / n)), flush=True)
print("  h = max sigma = %s (status %s, dual bound %s) (%.0fs)" % (h, r.status, hb, time.time() - t0), flush=True)
rp = solve_pair(); hp = -rp.fun if rp.x is not None else None
print("  h_pair = %s (status %s, dual bound %s) -> rigidity (R): %s (%.0fs)" % (hp, rp.status, getattr(rp, 'mip_dual_bound', None),
      (hp is not None and hp < h - 0.5) if rp.status == 0 else "undecided", time.time() - t0), flush=True)
if len(sys.argv) > 4: sys.exit(0)
rz = solve_z2(); z2 = -rz.fun if rz.x is not None else None
print("  z2 = max |S1|+|S2| = %s (status %s) ; 2h = %s -> box-zonotope certificate of (N): %s (%.0fs)" % (
      z2, rz.status, 2 * h, (z2 is not None and z2 >= 2 * h - 0.5), time.time() - t0), flush=True)
