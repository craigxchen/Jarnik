"""NUMERICAL: the quartic family Q4 (all 10 coordinate planes + hidden triples 134,135,145,234,235,245,345).
c_m = lam_m prod_{e != m}(t - r_me).  Gauge: lam5 = 1, r12 = 0, r34 = 1, t345 = 2.
Find members, check class and dimension; then estimate deg ev_{s0}."""
import numpy as np, itertools, sys, warnings
warnings.filterwarnings('ignore')
from levelsets import analyse
P = np.polynomial.polynomial
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
E = list(itertools.combinations(range(1, 6), 2))
HID = [(1,3,4), (1,3,5), (1,4,5), (2,3,4), (2,3,5), (2,4,5), (3,4,5)]
# unknown vector: r (8: all edges except (1,2),(3,4)), lam1..lam4, t_S for S != 345 (6)
FREE_E = [e for e in E if e not in [(1, 2), (3, 4)]]
HS = [S for S in HID if S != (3, 4, 5)]
def unpack(u):
    r = {(1, 2): 0.0, (3, 4): 1.0}
    for k, e in enumerate(FREE_E): r[e] = u[k]
    lam = list(u[8:12]) + [1.0]
    ts = {(3, 4, 5): 2.0}
    for k, S in enumerate(HS): ts[S] = u[12 + k]
    return r, lam, ts
def coords(r, lam):
    cs = []
    for m in range(1, 6):
        p = np.array([lam[m - 1]], complex)
        for e in range(1, 6):
            if e != m: p = P.polymul(p, [-r[tuple(sorted((m, e)))], 1])
        cs.append(p)
    return cs
def eqs(u):
    r, lam, ts = unpack(u); cs = coords(r, lam)
    out = []
    for S in HID:
        t = ts[S]; v = [P.polyval(t, cs[i - 1]) for i in S]
        out += [v[0] - v[1], v[1] - v[2]]
    return np.array(out)
def newton(u, fixed_idx, it=100):
    u = np.array(u, complex); var = [k for k in range(18) if k not in fixed_idx]
    for _ in range(it):
        f = eqs(u); J = np.zeros((14, len(var)), complex)
        for a, k in enumerate(var):
            du = np.zeros(18, complex); du[k] = 1e-7 * (1 + abs(u[k]))
            J[:, a] = (eqs(u + du) - eqs(u - du)) / (2 * du[k])
        st = np.linalg.lstsq(J, -f, rcond=None)[0]
        if not np.all(np.isfinite(st)): return None
        u[var] += st
        if np.linalg.norm(st) < 1e-13 * (1 + np.linalg.norm(u)): break
    return u if np.linalg.norm(eqs(u)) < 1e-10 else None
EXPECT = {}
for d, e in E: EXPECT[frozenset({d, e, 6})] = 1
for S in HID: EXPECT[frozenset(S)] = 1
EXPECT[frozenset({1, 2})] = 3
for p in [(1,3),(1,4),(1,5),(2,3),(2,4),(2,5)]: EXPECT[frozenset(p)] = 1
def genuine(u):
    r, lam, ts = unpack(u)
    cnt, bad, md, infb = analyse(coords(r, lam))
    return dict(cnt) == EXPECT and not bad and md > 1e-4 and not infb, cnt, bad, md
if __name__ == '__main__':
    found = 0
    for trial in range(int(sys.argv[2]) if len(sys.argv) > 2 else 200):
        u0 = rng.normal(size=18) + 1j * rng.normal(size=18)
        fixed = [0, 5, 8, 12]   # r13, r35?, lam1, t134 fixed (4 = dimension)
        u = newton(u0, fixed)
        if u is None: continue
        ok, cnt, bad, md = genuine(u)
        if ok:
            found += 1
            # Jacobian rank of the 14 equations w.r.t. all 18 unknowns
            J = np.zeros((14, 18), complex)
            for k in range(18):
                du = np.zeros(18, complex); du[k] = 1e-7 * (1 + abs(u[k]))
                J[:, k] = (eqs(u + du) - eqs(u - du)) / (2 * du[k])
            sv = np.linalg.svd(J, compute_uv=False)
            print('genuine member; Jacobian singular values min/max %.3g / %.3g, rank %d' % (sv[-1], sv[0], np.sum(sv > 1e-8 * sv[0])), flush=True)
            np.save('q4_member_%d.npy' % found, u)
            if found >= 3: break
        else:
            pass
    print('found', found)
