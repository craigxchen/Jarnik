"""NUMERICAL (referee): twisted cubics in P^4 meeting the 13 Kapranov planes of L7 Prop 3 through a random
general point, with a barrier unknown w:  w * Delta(u) = 1,  Delta = det(4x4 minor of [C0..C3]) *
prod_{i<j} (s_i - s_j) / scale  (all 13 contact parameters distinct; cubic spans a P^3).
Same gauge as cubic13_exist.py.  Accepted solutions are then tested for: meeting the 7 unmet planes,
the 15 lines L_ab, the 6 points p_a  (all must be avoided), and the predicted D_bc-intersections."""
import numpy as np, itertools, sys
exec(open('cubic13_exist.py').read().split('found = []')[0])
NU = 27
pairs = list(itertools.combinations(range(13), 2))
def Delta(u):
    C, s = unpack(u[:26])
    M = C.T                     # 5 x 4
    d = np.linalg.det(M[1:5, :])
    pr = np.prod([(s[i]-s[j]) for i, j in pairs])
    return d*pr
def G(u):
    return np.concatenate([F(u[:26]), [u[26]*Delta(u) - 1]])
def JG(u, h=1e-7):
    f0 = G(u); n = len(u); Jm = np.zeros((len(f0), n), complex)
    for k in range(n):
        d = np.zeros(n, complex); d[k] = h
        Jm[:, k] = (G(u + d) - G(u - d))/(2*h)
    return Jm
def point_test(C):
    best = np.inf
    for a in range(6):
        pa = pts[a]
        # c(s) proportional to pa:  all 2x2 minors of [c(s), pa] vanish; test at roots of one minor
        i0 = int(np.argmax(abs(pa)))
        for k in range(5):
            if k == i0: continue
            poly = np.array([C[3][k]*pa[i0] - C[3][i0]*pa[k], C[2][k]*pa[i0] - C[2][i0]*pa[k],
                             C[1][k]*pa[i0] - C[1][i0]*pa[k], C[0][k]*pa[i0] - C[0][i0]*pa[k]])
            for r in np.roots(poly):
                x = curve(C, r); x = x/np.linalg.norm(x)
                y = pa/np.linalg.norm(pa)
                best = min(best, np.linalg.norm(x - (np.vdot(y, x))*y))
            break
    return best
sols = []
for trial in range(300):
    u = np.concatenate([rng.normal(size=26) + 1j*rng.normal(size=26), [0.0]])
    u[26] = 1/Delta(u)
    ok = False
    for it in range(100):
        f = G(u); nf = np.linalg.norm(f)
        if nf < 1e-11: ok = True; break
        du = np.linalg.lstsq(JG(u), -f, rcond=None)[0]
        lam = 1.0
        while lam > 1e-5 and np.linalg.norm(G(u + lam*du)) > nf: lam /= 2
        u = u + lam*du
        if not np.all(np.isfinite(u)) or np.linalg.norm(u[:26]) > 1e6: break
    if not ok: continue
    C, s = unpack(u[:26])
    sv = np.linalg.svd(J(u[:26]), compute_uv=False)
    rC = np.linalg.svd(C, compute_uv=False)
    sd = min(abs(a - b) for a, b in itertools.combinations(s, 2))
    op = min(meets(C, ANN[p]) for p in UNMET); ln = min(meets(C, LANN[l]) for l in LANN)
    pt = point_test(C)
    sols.append(np.concatenate([q, u]))
    print('trial %3d |F| %.1e  sv_min/sv_max(J_F) %.2e  C-indep %.2e  min|si-sj| %.2e  other-planes %.2e  lines %.2e  points %.2e'
          % (trial, np.linalg.norm(F(u[:26])), sv[-1]/sv[0], rC[-1]/rC[0], sd, op, ln, pt), flush=True)
    if len(sols) >= 5: break
print('accepted', len(sols), 'of', trial + 1)
np.save('cubic13_sols_%s.npy' % (sys.argv[1] if len(sys.argv) > 1 else '0'), np.array(sols))
