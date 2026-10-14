"""EVIDENCE ONLY (floating point).  Paley-12, T = 1, L = 2, m = 8 (zonotope-free, LP-subcritical).
Every 8-point section Y (base + 7 of the 11 coordinates) has dim V(A_Y, 8) = 1 and its octic g_Y is irreducible and
indefinite (section_forms.py, section_sign.py), so each section alone has nondegenerate real zeros.
By Prop 3.3, TC(A,8) is contained in the common zero set Z of all 330 section octics g_Y(u_Y).
Here we search numerically for real points of Z on the unit sphere of R^11 and measure how degenerate they are.
"""
import itertools, sys, time
import numpy as np
import flint
from scipy.optimize import least_squares
from tclib import box_points, cube_monomials, mon_eval

q = 11; D = [0, 1, 3, 4, 5, 9]; L = 2
cols = [sorted((i + x) % q for x in D) for i in range(q)]


def kernel(pts, m, d):
    mons = [a for a in cube_monomials(d) if sum(a) < m]
    M = flint.fmpz_mat(len(mons), len(pts), [mon_eval(k, a) for a in mons for k in pts])
    X, nullity = M.nullspace()
    return [[int(X[i, j]) for i in range(len(pts))] for j in range(nullity)]


t0 = time.time()
cache = {}
sections = []
for Y in itertools.combinations(range(q), 7):
    can, t = min((tuple(sorted((y + s) % q for y in Y)), s) for s in range(q))
    if can not in cache:
        idx = {y: i for i, y in enumerate(can)}
        rc = [[idx[x] for x in c if x in idx] for c in cols]; rc = [c for c in rc if c]
        pts = box_points(7, [(c, L) for c in rc])
        c = np.array(kernel(pts, 8, 7)[0], dtype=float)
        c /= np.abs(c).sum()
        cache[can] = (np.array(pts, dtype=float), c)
    K, c = cache[can]
    # coordinates: section coordinate i corresponds to point can[i] = (Y_j + t) -> original Y_j = can[i] - t
    orig = [(y - t) % q for y in can]
    sections.append((orig, K, c))
print('sections %d, classes %d (%.0fs)' % (len(sections), len(cache), time.time() - t0), flush=True)


def resid(v):
    u = v / np.linalg.norm(v)
    out = []
    for orig, K, c in sections:
        out.append(((K @ u[orig]) ** 8) @ c)
    return np.array(out)


def degeneracy(u):
    # smallest |<v,u>| over v in {e_x, e_x - e_y, e_x + e_y - e_z, e_x - e_y - e_z + e_w ...} small set, relative
    vals = [abs(x) for x in u]
    for x, y in itertools.combinations(range(q), 2):
        vals.append(abs(u[x] - u[y])); vals.append(abs(u[x] + u[y]))
    for x, y, z in itertools.permutations(range(q), 3):
        if x < y:
            vals.append(abs(u[x] + u[y] - u[z]))
    return min(vals)


rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
ntr = int(sys.argv[2]) if len(sys.argv) > 2 else 200
results = []
for tr in range(ntr):
    v0 = rng.standard_normal(q)
    sol = least_squares(resid, v0, xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=4000)
    u = sol.x / np.linalg.norm(sol.x)
    r = np.linalg.norm(resid(u))
    dg = degeneracy(u)
    results.append((r, dg))
    print('trial %d  |resid| = %.3e  degeneracy = %.3e  u = %s' % (tr, r, dg, np.array2string(u, precision=4, max_line_width=200)), flush=True)
results.sort()
print('smallest residuals with their degeneracy:', results[:10])
