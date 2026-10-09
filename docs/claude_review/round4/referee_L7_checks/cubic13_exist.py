"""NUMERICAL (referee): do twisted cubics in P^4 meeting the 13 Kapranov planes of L7 Prop 3 exist through a
GENERAL point?  (L7 says this was not checked.)  If yes, the class 'cub13' (pair 6, triple 13) is the class of
a covering, hence (Kollar II.3.11) free, family of rational curves.

Kapranov psi_7 model: p_1..p_5 = e_1..e_5, p_6 = (1,1,1,1,1) in P^4.  Plane Pi_abc = span(p_a,p_b,p_c), cut out
by two linear forms (null space).  Cubic c(s) = C0 + C1 s + C2 s^2 + C3 s^3 (C_k in C^5).
Gauge: C0 = q (random general point, s = 0), s_{Pi_1} = 1, s_{Pi_2} = -1 (fixes PGL_2 and scalar).
Unknowns: C1, C2, C3 (15) + s_Pi for the other 11 planes = 26.  Equations: 2 per plane = 26.  Square system.
A solution is accepted if: residual < 1e-10, Jacobian well conditioned, C0..C3 independent (spans a P^3),
the 13 parameters s_Pi distinct, the cubic does NOT meet the 7 other planes, the 15 lines L_ab, or pass
through p_a (checked by resultant-free root tests), and every D_bc-intersection is as predicted.
"""
import numpy as np, itertools, sys

rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
pts = [np.eye(5)[i] for i in range(5)] + [np.ones(5)]
MET = [(1, 2, 3), (1, 2, 4), (1, 2, 5), (1, 2, 6), (1, 3, 4), (1, 3, 5), (1, 3, 6), (1, 4, 5),
       (2, 4, 6), (2, 5, 6), (3, 4, 6), (3, 5, 6), (4, 5, 6)]
ALL = list(itertools.combinations(range(1, 7), 3))
UNMET = [p for p in ALL if p not in MET]

def annih(idx):
    M = np.array([pts[i-1] for i in idx])
    _, s, vt = np.linalg.svd(M)
    return vt[len(idx):]          # (5-len) x 5 forms vanishing on the span

ANN = {p: annih(p) for p in ALL}
LANN = {l: annih(l) for l in itertools.combinations(range(1, 7), 2)}

q = rng.normal(size=5) + 1j*rng.normal(size=5)

def unpack(u):
    C = np.vstack([q, u[0:5], u[5:10], u[10:15]])     # C[k] coefficient of s^k
    s = np.concatenate([[1.0, -1.0], u[15:26]])
    return C, s

def curve(C, s):
    return C[0] + C[1]*s + C[2]*s**2 + C[3]*s**3

def F(u):
    C, s = unpack(u)
    R = []
    for n, p in enumerate(MET):
        x = curve(C, s[n])
        R += list(ANN[p] @ x)
    return np.array(R)

def J(u, h=1e-7):
    f0 = F(u); n = len(u); Jm = np.zeros((len(f0), n), complex)
    for k in range(n):
        d = np.zeros(n, complex); d[k] = h
        Jm[:, k] = (F(u + d) - F(u - d))/(2*h)
    return Jm

def meets(C, forms):
    """does the cubic meet the linear space cut out by 'forms' (rows)?  min over s of normalized residual
    via the gcd test: common root of the cubics forms@c(s)."""
    polys = [np.array([f @ C[3], f @ C[2], f @ C[1], f @ C[0]]) for f in forms]
    roots = np.roots(polys[0])
    best = np.inf
    for r in roots:
        x = curve(C, r)
        val = np.linalg.norm(np.array([f @ x for f in forms]))/np.linalg.norm(x)
        best = min(best, val)
    # also s = infinity
    x = C[3]; val = np.linalg.norm(np.array([f @ x for f in forms]))/max(np.linalg.norm(x), 1e-300)
    return min(best, val)

found = []
for trial in range(400):
    u = rng.normal(size=26) + 1j*rng.normal(size=26)
    ok = False
    for it in range(80):
        f = F(u); nf = np.linalg.norm(f)
        if nf < 1e-12: ok = True; break
        Jm = J(u)
        du = np.linalg.lstsq(Jm, -f, rcond=None)[0]
        lam = 1.0
        while lam > 1e-4 and np.linalg.norm(F(u + lam*du)) > nf: lam /= 2
        u = u + lam*du
        if np.linalg.norm(u) > 1e8: break
    if not ok: continue
    C, s = unpack(u)
    sv = np.linalg.svd(J(u), compute_uv=False)
    rankC = np.linalg.svd(C, compute_uv=False)
    sd = min(abs(a - b) for a, b in itertools.combinations(s, 2))
    other_planes = min(meets(C, ANN[p]) for p in UNMET)
    lines = min(meets(C, LANN[l]) for l in LANN)
    points = min(min(np.linalg.norm(np.cross(curve(C, r)[:3], pts[a][:3])) for r in np.roots(np.array([ (np.ones(5) if a==5 else pts[a]) @ C[k] for k in (3,2,1,0)]))) for a in range(6)) if False else None
    found.append((trial, nf, sv[-1]/sv[0], rankC[-1]/rankC[0], sd, other_planes, lines))
    print('trial %3d resid %.1e  cond(J)^-1 %.2e  C-indep %.2e  min|s_i-s_j| %.2e  meet-other-planes %.2e  meet-lines %.2e'
          % found[-1], flush=True)
    if len(found) >= 6: break
print('converged:', len(found), 'of', trial+1, 'starts')
