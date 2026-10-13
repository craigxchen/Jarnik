"""NUMERICAL: Q4 family with gauge t345 = inf (lam3=lam4=lam5=1), r34 = 0, r35 = 1.
Unknowns u: r12 r13 r14 r15 r23 r24 r25 r45 lam1 lam2 t134 t234 t135 t235 t145 t245 (16); 12 equations."""
import numpy as np, sys, warnings, itertools
warnings.filterwarnings('ignore')
from levelsets import analyse
P = np.polynomial.polynomial
NAMES = 'r12 r13 r14 r15 r23 r24 r25 r45 lam1 lam2 t134 t234 t135 t235 t145 t245'.split()
def coords(u):
    d = dict(zip(NAMES, u)); r = {(3, 4): 0.0, (3, 5): 1.0}
    for n in NAMES[:8]: r[(int(n[1]), int(n[2]))] = d[n]
    lam = [d['lam1'], d['lam2'], 1, 1, 1]
    cs = []
    for m in range(1, 6):
        p = np.array([lam[m - 1]], complex)
        for e in range(1, 6):
            if e != m: p = P.polymul(p, [-r[tuple(sorted((m, e)))], 1])
        cs.append(p)
    return cs, d
CONDS = [(3, 4, 't134'), (3, 4, 't234'), (3, 5, 't135'), (3, 5, 't235'), (4, 5, 't145'), (4, 5, 't245'),
         (1, 3, 't134'), (1, 3, 't135'), (1, 4, 't145'), (2, 3, 't234'), (2, 3, 't235'), (2, 4, 't245')]
def eqs(u):
    cs, d = coords(u)
    return np.array([P.polyval(d[t], cs[i - 1]) - P.polyval(d[t], cs[j - 1]) for i, j, t in CONDS])
def newton(u, fixed, it=200):
    u = np.array(u, complex); var = [k for k in range(16) if k not in fixed]
    for k_ in range(it):
        f = eqs(u); J = np.zeros((12, 12), complex)
        for a, k in enumerate(var):
            du = np.zeros(16, complex); du[k] = 1e-7 * (1 + abs(u[k]))
            J[:, a] = (eqs(u + du) - eqs(u - du)) / (2 * du[k])
        try: st = np.linalg.solve(J, -f)
        except np.linalg.LinAlgError: return None
        if not np.all(np.isfinite(st)): return None
        nst = np.linalg.norm(st)
        if nst > 1.0: st *= 1.0 / nst   # damping
        u[var] += st
        if nst < 1e-13 * (1 + np.linalg.norm(u)): break
    return u if np.linalg.norm(eqs(u)) < 1e-10 * (1 + np.linalg.norm(u)) ** 4 else None
E = list(itertools.combinations(range(1, 6), 2))
HID = [(1,3,4), (1,3,5), (1,4,5), (2,3,4), (2,3,5), (2,4,5), (3,4,5)]
EXPECT = {}
for dd, ee in E: EXPECT[frozenset({dd, ee, 6})] = 1
for S in HID: EXPECT[frozenset(S)] = 1
EXPECT[frozenset({1, 2})] = 3
for p in [(1,3),(1,4),(1,5),(2,3),(2,4),(2,5)]: EXPECT[frozenset(p)] = 1
def genuine(u):
    cs, d = coords(u)
    cnt, bad, md, infb = analyse(cs)
    return dict(cnt) == EXPECT and not bad and md > 1e-4 and not infb, cnt, bad, md
if __name__ == '__main__':
    rng = np.random.default_rng(int(sys.argv[1]))
    nconv = 0; ngen = 0
    for trial in range(int(sys.argv[2])):
        u0 = rng.normal(size=16) * 1.5 + 1j * rng.normal(size=16) * 1.5
        fixed = list(rng.choice(16, size=4, replace=False))
        u = newton(u0, fixed)
        if u is None: continue
        nconv += 1
        ok, cnt, bad, md = genuine(u)
        if ok:
            ngen += 1; np.save('q4b_member_%s_%d.npy' % (sys.argv[1], ngen), u)
            print('GENUINE', np.round(u, 4), flush=True)
            if ngen >= 2: break
    print('converged', nconv, 'genuine', ngen)
