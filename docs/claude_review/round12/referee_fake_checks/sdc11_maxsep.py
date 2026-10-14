"""Evidence: maximise the separation over the real common zero set of the 5 sdc11 octics (T=1, L=2, m=8, base inf).
SLSQP with equality constraints F(u)=0 (normalised), u_0 = 1, |u_x| <= 1, pair gaps >= t; maximise t.
Starts: random, and the two certified zeros."""
import json, sys, time, numpy as np
from scipy.optimize import minimize
B = json.load(open('sdc11_inv_basis.json'))
K = np.array(B['K'], dtype=float); Cs = [np.array(c, dtype=float) for c in B['C']]
d = 11; m = 8
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
U = rng.uniform(-1, 1, (200, d)); U[:, 0] = 1
L = U @ K.T
scale = [np.sqrt(np.mean((L**m @ c)**2)) for c in Cs]
Cn = np.array([c/s for c, s in zip(Cs, scale)])          # 5 x |A|
def F(u): return Cn @ (K @ u)**m
def JF(u): return (Cn * (m*(K @ u)**(m-1))) @ K
pairs = [(i, j) for i in range(d+1) for j in range(i+1, d+1)]
def z(u): return np.concatenate([[0.0], u])
def solve(u0, order):
    # variables x = (u_1..u_10, t) with u_0 = 1 fixed; ordering constraint z_{order[a+1]} - z_{order[a]} >= t
    def unpack(x): return np.concatenate([[1.0], x[:d-1]]), x[-1]
    cons = [{'type': 'eq', 'fun': lambda x: F(unpack(x)[0]), 'jac': lambda x: np.hstack([JF(unpack(x)[0])[:, 1:], np.zeros((5, 1))])}]
    def gfun(x):
        u, t = unpack(x); zz = z(u)
        return np.array([zz[order[a+1]] - zz[order[a]] - t for a in range(d)])
    def gjac(x):
        J = np.zeros((d, d)); 
        for a in range(d):
            for idx, sgn in ((order[a+1], 1), (order[a], -1)):
                if idx >= 2: J[a, idx-2] += sgn       # z index idx -> u index idx-1 -> x index idx-2
            J[a, -1] = -1
        return J
    cons.append({'type': 'ineq', 'fun': gfun, 'jac': gjac})
    x0 = np.concatenate([u0[1:], [0.0]])
    res = minimize(lambda x: -x[-1], x0, jac=lambda x: np.concatenate([np.zeros(d-1), [-1.0]]), constraints=cons,
                   bounds=[(-1, 1)]*(d-1) + [(-1, 1)], method='SLSQP', options={'maxiter': 500, 'ftol': 1e-14})
    u, t = unpack(res.x)
    zz = np.sort(z(u)); sep = np.min(np.diff(zz))
    return u, sep, np.max(np.abs(F(u)))
best = []
starts = [np.array(x) for x in np.load('../fake_checks/sdc11_landzeros.npy')]
for _ in range(150):
    u = rng.uniform(-1, 1, d); u[0] = 1; starts.append(u)
t0 = time.time()
for i, u0 in enumerate(starts):
    order = list(np.argsort(z(u0)))
    try:
        u, sep, res = solve(u0, order)
    except Exception as e:
        continue
    if res < 1e-10:
        best.append((sep, u.tolist(), res))
        print('start %2d: feasible zero, sep = %.5f, max|F| = %.1e (%.0fs)' % (i, sep, res, time.time()-t0), flush=True)
best.sort(key=lambda x: -x[0])
json.dump(best[:10], open('sdc11_maxsep_%s.json' % (sys.argv[1] if len(sys.argv) > 1 else '0'), 'w'))
print('best separations found:', [round(b[0], 5) for b in best[:10]])
