"""NUMERICAL: Q4 via the projection structure (cf. cub13).  Gauge t345 = inf (lam3=lam4=lam5=1),
r34 = 0, r35 = 1.  Free: t234, t235, t245, r23, r24, r25, r45 (7).
(A),(B),(C): roots t234, t235, t245 of B34, B35, B45  -> linear 3x3 system for (r13, r14, r15).
(D): lam2 (t - r12) P_k(t) = R_k(t) at t234, t235 (P=(t-r24)(t-r25), R=(t-r13)(t)(t-1)) and
     at t245 (P=(t-r23)(t-r25), R=(t-r14) t (t-r45)): 3 eqs linear in (lam2, mu=lam2*r12) -> 1 condition.
Lift (label 1): t134, t135, t145 = other roots of B34, B35, B45;
     lam1 (t-r12)(t-r14)(t-r15) = (t-r23) t (t-1) at t134, t135;  lam1 (t-r12)(t-r13)(t-r15) = (t-r24) t (t-r45) at t145.
Unknowns solved by Newton: 3 of the 7 free parameters; equations: (D)-consistency + 2 lift conditions."""
import numpy as np, sys, warnings, itertools
warnings.filterwarnings('ignore')
from levelsets import analyse
P = np.polynomial.polynomial
FREE = ['t234', 't235', 't245', 'r23', 'r24', 'r25', 'r45']
def lin(r): return np.array([-r, 1], complex)
def pm(*ps):
    out = np.array([1], complex)
    for p in ps: out = P.polymul(out, p)
    return out
def solve_all(f):
    t234, t235, t245, r23, r24, r25, r45 = [f[k] for k in FREE]
    r34, r35 = 0.0, 1.0
    # B34(t) = (t-r13)(t-r23)(t-r35) - (t-r14)(t-r24)(t-r45) ; linear in r13, r14:
    # (t-r13)A - (t-r14)Bq = t A - t Bq - r13 A + r14 Bq
    def rowB(t, A, Bq, which):   # coefficients for (r13, r14, r15) and rhs
        pass
    A34 = lambda t: (t - r23) * (t - r35); B34q = lambda t: (t - r24) * (t - r45)
    A35 = lambda t: (t - r23) * (t - r34); B35q = lambda t: (t - r25) * (t - r45)
    A45 = lambda t: (t - r24) * (t - r34); B45q = lambda t: (t - r25) * (t - r35)
    M = np.array([[-A34(t234), B34q(t234), 0], [-A35(t235), 0, B35q(t235)], [0, -A45(t245), B45q(t245)]], complex)
    rhs = -np.array([t234 * (A34(t234) - B34q(t234)), t235 * (A35(t235) - B35q(t235)), t245 * (A45(t245) - B45q(t245))])
    r13, r14, r15 = np.linalg.solve(M, rhs)
    # (D): lam2 (t - r12) Pk = Rk  ->  lam2 t Pk - mu Pk = Rk
    rows = []; rr = []
    for t, Pk, Rk in [(t234, (t234 - r24) * (t234 - r25), (t234 - r13) * t234 * (t234 - 1)),
                      (t235, (t235 - r24) * (t235 - r25), (t235 - r13) * t235 * (t235 - 1)),
                      (t245, (t245 - r23) * (t245 - r25), (t245 - r14) * t245 * (t245 - r45))]:
        rows.append([t * Pk, -Pk]); rr.append(Rk)
    rows = np.array(rows); rr = np.array(rr)
    lam2, mu = np.linalg.solve(rows[:2], rr[:2])
    eqD = rows[2] @ np.array([lam2, mu]) - rr[2]
    r12 = mu / lam2
    # other roots
    def quad(Aq, Bq_, t0, r1x, r1y):
        # B(t) = (t - r1x) A(t) - (t - r1y) Bq(t): build polynomial
        pass
    B34 = P.polysub(pm(lin(r13), lin(r23), lin(r35)), pm(lin(r14), lin(r24), lin(r45)))
    B35 = P.polysub(pm(lin(r13), lin(r23), lin(r34)), pm(lin(r15), lin(r25), lin(r45)))
    B45 = P.polysub(pm(lin(r14), lin(r24), lin(r34)), pm(lin(r15), lin(r25), lin(r35)))
    def other(B, t0):
        B = np.trim_zeros(B, 'b') if np.abs(B[-1]) < 1e-14 else B
        B = B[:3]   # quadratic (cubic terms cancel)
        return -B[1] / B[2] - t0
    t134, t135, t145 = other(B34, t234), other(B35, t235), other(B45, t245)
    l1a = ((t134 - r23) * t134 * (t134 - 1)) / ((t134 - r12) * (t134 - r14) * (t134 - r15))
    l1b = ((t135 - r23) * t135 * (t135 - 1)) / ((t135 - r12) * (t135 - r14) * (t135 - r15))
    l1c = ((t145 - r24) * t145 * (t145 - r45)) / ((t145 - r12) * (t145 - r13) * (t145 - r15))
    r = {(1,2): r12, (1,3): r13, (1,4): r14, (1,5): r15, (2,3): r23, (2,4): r24, (2,5): r25, (3,4): r34, (3,5): r35, (4,5): r45}
    lam = [l1a, lam2, 1, 1, 1]
    return np.array([eqD, l1a - l1b, l1a - l1c]), r, lam
def coords(r, lam):
    cs = []
    for m in range(1, 6):
        p = np.array([lam[m - 1]], complex)
        for e in range(1, 6):
            if e != m: p = P.polymul(p, lin(r[tuple(sorted((m, e)))]))
        cs.append(p)
    return cs
HID = [(1,3,4), (1,3,5), (1,4,5), (2,3,4), (2,3,5), (2,4,5), (3,4,5)]
EXPECT = {}
for d_, e_ in itertools.combinations(range(1, 6), 2): EXPECT[frozenset({d_, e_, 6})] = 1
for S in HID: EXPECT[frozenset(S)] = 1
EXPECT[frozenset({1, 2})] = 3
for p_ in [(1,3),(1,4),(1,5),(2,3),(2,4),(2,5)]: EXPECT[frozenset(p_)] = 1
def F(x, base, unk):
    f = dict(base); f.update(dict(zip(unk, x)))
    return solve_all(f)[0]
def newton(x, base, unk):
    x = np.array(x, complex)
    for _ in range(100):
        fx = F(x, base, unk); J = np.zeros((3, 3), complex)
        for k in range(3):
            dx = np.zeros(3, complex); dx[k] = 1e-7 * (1 + abs(x[k]))
            J[:, k] = (F(x + dx, base, unk) - F(x - dx, base, unk)) / (2 * dx[k])
        st = np.linalg.solve(J, -fx)
        if not np.all(np.isfinite(st)): return None
        n = np.linalg.norm(st)
        if n > 1: st /= n
        x = x + st
        if n < 1e-13 * (1 + np.linalg.norm(x)): break
    return x if np.linalg.norm(F(x, base, unk)) < 1e-10 else None
if __name__ == '__main__':
    rng = np.random.default_rng(int(sys.argv[1])); found = 0; conv = 0
    for trial in range(int(sys.argv[2])):
        base = {k: complex(rng.normal() * 2, rng.normal() * 2) for k in FREE}
        unk = list(rng.choice(FREE, size=3, replace=False))
        try:
            x = newton([base[k] for k in unk], base, unk)
        except Exception:
            x = None
        if x is None: continue
        conv += 1
        f = dict(base); f.update(dict(zip(unk, x)))
        try:
            _, r, lam = solve_all(f)
            cnt, bad, md, infb = analyse(coords(r, lam))
        except Exception:
            continue
        if dict(cnt) == EXPECT and not bad and md > 1e-4 and not infb:
            found += 1; np.save('q4c_member_%s_%d.npy' % (sys.argv[1], found), np.array([f[k] for k in FREE]))
            print('GENUINE Q4 member; min contact distance %.3g' % md, {k: np.round(v, 4) for k, v in f.items()}, flush=True)
            if found >= 3: break
    print('converged', conv, 'genuine', found)
