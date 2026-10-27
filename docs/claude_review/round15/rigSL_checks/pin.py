# pin.py NAME K [rho] [maxdepth] -- angle pinning by unit-D monomials (round 15, Rig^SL).
# Claim checked: every phi in [0,pi/4]^r (fundamental domain of the symmetry group, Lemma 2.1 of rigSL.md) with
#     |sin(<v,phi> - s pi/4)| >= nu_s e^(-w(v)/2)    for all 1 <= ||v||_1 <= K, all s          (*)
# (the label-free part of (1.2): m_(v,s) >= 1) equals the folded genuine angle vector g.
# (a) local: T = {(v,s): |D_s(Z_v)| = 1} (constraints of (*) tight at g); feasible directions satisfy sigma<v,d> >= 0;
#     exact certificate that this cone is {0}: rational lambda > 0 on a support U, sum lambda_k sigma_k v_k = 0, rank U = r.
# (b) global: box branch-and-bound excludes [0,pi/4]^r minus the cube |phi - g|_inf <= rho (float, margin 1e-9 in
#     favour of feasibility: exclusion is sound).
# (c) inside the cube every tight constraint is a half-space: K*rho < 2*alpha_min (alpha = arcsin(nu_s e^(-w/2))).
import sys, math
from fractions import Fraction
import numpy as np
from scipy.optimize import linprog
import sympy
from sl import *

def fold(x):
    x = x % (math.pi / 2)
    return min(x, math.pi / 2 - x), (1 if x <= math.pi / 4 else -1)

def run(name, K, rho=1e-5, maxdepth=200, verbose=True):
    S = setup(name, K); r = S.r
    g = []; sg = []
    for (a, b) in S.pis:
        x, e = fold(math.atan2(b, a)); g.append(x); sg.append(e)
    # folded system: column j conjugated iff sg[j] = -1; then D_s of the folded system = D_s of eps*v (up to sign)
    # folded Gaussian primes pi'_j (associate / conjugate of pi_j) with arg pi'_j = g_j mod pi
    pisf = []
    for j, (a, b) in enumerate(S.pis):
        cands = []
        for z0 in ((a, b), (a, -b)):
            z = z0
            for _ in range(4):
                cands.append(z); z = (-z[1], z[0])
        z = [c for c in cands if abs((math.atan2(c[1], c[0]) - g[j] + math.pi / 2) % math.pi - math.pi / 2) < 1e-12]
        assert z; pisf.append(z[0])
    T = []
    for i, v in enumerate(S.V):
        X, Y = Zv(pisf, v); D = (Y, X - Y, X, X + Y)
        for s in range(4):
            N = S.N[i]; alpha = math.asin((2 ** -0.5 if s % 2 else 1.0) / math.sqrt(N))
            th = sum(a * b for a, b in zip(v, g)) - s * math.pi / 4
            k = round(th / math.pi); dev = th - k * math.pi
            assert (abs(D[s]) == 1) == (abs(abs(dev) - alpha) < 1e-9), (v, s, D[s], dev, alpha)
            if abs(D[s]) == 1:
                sigma = 1 if dev > 0 else -1
                T.append((tuple(sigma * t for t in v), alpha, v, s))
    out0 = pisf
    W = np.array([w for w, _, _, _ in T], dtype=float)
    out = dict(tight=len(T))
    # (a) LP: lambda in [0,1], W^T lambda = 0, maximise sum lambda
    res = linprog(-np.ones(len(T)), A_eq=W.T, b_eq=np.zeros(r), bounds=[(0, 1)] * len(T), method="highs")
    lam = res.x; U = [k for k in range(len(T)) if lam[k] > 1e-7]
    WU = sympy.Matrix([list(T[k][0]) for k in U]).T          # r x |U|
    rankU = WU.rank(); ns = WU.nullspace()
    cert = None
    if rankU == r and ns:
        B = sympy.Matrix.hstack(*ns)                          # |U| x d
        coef, *_ = np.linalg.lstsq(np.array(B, dtype=float), lam[U], rcond=None)
        for den in (1, 10, 100, 1000, 10 ** 4, 10 ** 6):
            c = [sympy.Rational(Fraction(x).limit_denominator(den)) for x in coef]
            l = B * sympy.Matrix(c)
            if all(x > 0 for x in l) and WU * l == sympy.zeros(r, 1):
                cert = [l[i] for i in range(len(U))]; break
    out["local_cone_trivial"] = cert is not None
    out["cert"] = (cert, [T[k][2:] for k in U]) if cert else None
    amin = min(a for _, a, _, _ in T)
    out["half_space_ok"] = K * rho < 2 * amin
    # (b) global exclusion
    boxes = [([0.0] * r, [math.pi / 4] * r)]; depth = 0
    def inside(lo, hi): return all(g[j] - rho <= lo[j] and hi[j] <= g[j] + rho for j in range(r))
    done_in = []
    while boxes and depth < maxdepth and len(boxes) < 200000:
        nb = []
        for lo, hi in boxes:
            if inside(lo, hi): done_in.append((lo, hi)); continue
            j = max(range(r), key=lambda k: hi[k] - lo[k]); mid = (lo[j] + hi[j]) / 2
            for a, b in ((lo[j], mid), (mid, hi[j])):
                l2 = list(lo); h2 = list(hi); l2[j] = a; h2[j] = b
                if min(budgets_box(S, l2, h2)) >= 0: nb.append((l2, h2))
        boxes = nb; depth += 1
    out["global_excluded"] = not boxes
    out["boxes_left"] = len(boxes); out["boxes_inside"] = len(done_in); out["depth"] = depth
    out["g"] = g; out["alpha_min"] = amin
    if verbose: print(name, K, {k: v for k, v in out.items() if k != "cert"}, flush=True)
    return out

if __name__ == '__main__':
    name, K = sys.argv[1], int(sys.argv[2])
    rho = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-5
    o = run(name, K, rho)
    if o["cert"]:
        cert, sup = o["cert"]
        print("  Farkas certificate (lambda, (v, s)):")
        for l, (v, s) in zip(cert, sup): print("   ", l, v, s)
    print("PINNED" if o["local_cone_trivial"] and o["global_excluded"] and o["half_space_ok"] else "NOT PROVED")
