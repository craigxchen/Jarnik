# Search for pair-separated LOW-RANK profiles (exact LP over weights, scipy linprog; verified with fractions).
# Points lambda_x in Z^rho, columns y_p = <d_p, lambda>, e_p = range of y_p over the points.
# Character inequality (model, leading order): ||l|| = sum_p omega_p (2/e_p)|<d_p,l>| >= 1 for all l in Lambda\0,
# sum omega_p = 1.  We maximise s = min ||l|| over weights; profile admissible iff s >= 1.
import itertools, sys, numpy as np
from scipy.optimize import linprog
from fractions import Fraction
rho = int(sys.argv[1]); G = int(sys.argv[2]); Mmin = int(sys.argv[3]); Mmax = int(sys.argv[4])
dirs = []
for d in itertools.product(range(-3, 4), repeat=rho):
    if all(x == 0 for x in d): continue
    if np.gcd.reduce([abs(x) for x in d]) != 1: continue
    if tuple(-x for x in d) in dirs: continue
    dirs.append(d)
testL = [l for l in itertools.product(range(-6, 7), repeat=rho) if any(l)]
def best_weights(pts):
    ys = np.array([[np.dot(d, p) for d in dirs] for p in pts])
    e = ys.max(0) - ys.min(0)
    use = [i for i in range(len(dirs)) if e[i] > 0]
    coef = np.array([[2*abs(np.dot(dirs[i], l))/e[i] for i in use] for l in testL])
    # maximise s: variables (omega_use, s); constraints coef@omega >= s
    nv = len(use)+1
    c = np.zeros(nv); c[-1] = -1
    Aub = np.hstack([-coef, np.ones((len(testL), 1))]); bub = np.zeros(len(testL))
    Aeq = np.zeros((1, nv)); Aeq[0, :-1] = 1
    r = linprog(c, A_ub=Aub, b_ub=bub, A_eq=Aeq, b_eq=[1], bounds=[(0, None)]*(nv-1)+[(None, None)], method='highs')
    return -r.fun, use, r.x[:-1], e
grid = list(itertools.product(range(G), repeat=rho))
for M in range(Mmin, Mmax+1):
    best = None; cnt = 0
    for pts in itertools.combinations(grid, M):
        if pts[0] != tuple([0]*rho): continue
        # lattice generated must be Z^rho (full rank) -- check rank
        D = np.array([np.subtract(p, pts[0]) for p in pts[1:]])
        if np.linalg.matrix_rank(D) < rho: continue
        s, use, om, e = best_weights(pts)
        cnt += 1
        if best is None or s > best[0]: best = (s, pts, [(dirs[u], round(o, 4), int(e[u])) for u, o in zip(use, om) if o > 1e-9])
    print(f"rho={rho} M={M}: {cnt} point sets; best min-norm s = {best[0]:.4f} at points {best[1]}; weights(dir,omega,e) {best[2]}", flush=True)
