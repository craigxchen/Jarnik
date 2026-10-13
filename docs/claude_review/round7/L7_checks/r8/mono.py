# Monodromy computation of deg(ev) for a Mobius type: number of members through a general point.
import numpy as np, sys, itertools
from fwd import to_mobius
import mob
def member_to_X(V, E, t0):
    n = 7; m = len(E)
    s = V[:m]; z = V[m:2*m]; C = V[2*m:].reshape(n, 4)
    g = [to_mobius(C[a]) for a in range(n)]
    y = np.array([(g[a][0,0]*t0 + g[a][0,1])/(g[a][1,0]*t0 + g[a][1,1]) for a in range(n)])
    p = np.zeros(n, complex); q = np.zeros(n, complex)
    for a in range(n):
        A, B, Cc, D = g[a][0,0], g[a][0,1], g[a][1,0], g[a][1,1]
        p[a] = -Cc/(D + t0*Cc); q[a] = A/(D + t0*Cc)
    lam = 1/(q[0] + p[0]*y[0]); mu = -lam*p[0]
    p2 = lam*p + mu; q2 = lam*q - mu*y
    zz = np.array([(g[E[i][0]][0,0]*s[i] + g[E[i][0]][0,1])/(g[E[i][0]][1,0]*s[i] + g[E[i][0]][1,1]) for i in range(m)])
    X = np.concatenate([p2[1:], q2[1:], zz])
    return X, y
def track(X, yA, yB, E, h0=0.02):
    s = 0.0; h = h0; dy = yB - yA
    while s < 1.0:
        h = min(h, 1.0 - s)
        # RK4 predictor
        def vel(X, s):
            JX, Jy = mob.J(X, yA + s*dy, E, wrt_y=True)
            return -np.linalg.solve(JX, Jy @ dy)
        try:
            k1 = vel(X, s); k2 = vel(X + h/2*k1, s + h/2); k3 = vel(X + h/2*k2, s + h/2); k4 = vel(X + h*k3, s + h)
        except np.linalg.LinAlgError:
            return None
        Xp = X + h/6*(k1 + 2*k2 + 2*k3 + k4)
        y1 = yA + (s + h)*dy
        ok = False
        for it in range(6):
            f = mob.F(Xp, y1, E)
            dX = np.linalg.solve(mob.J(Xp, y1, E), -f)
            Xp = Xp + dX
            if np.linalg.norm(dX) < 1e-11*(1 + np.linalg.norm(Xp)): ok = True; break
            if it == 0 and np.linalg.norm(dX) > 1e-3*(1+np.linalg.norm(Xp)): break
        if ok and it <= 3:
            X = Xp; s += h; h = min(h*1.5, 0.1)
        else:
            h /= 2
            if h < 1e-9: return None
    return X
def monodromy(X0, y0, E, rng, nloops=40, scale=1.0, verbose=True):
    sols = [X0]; stale = 0
    for L in range(nloops):
        y1 = y0 + scale*(rng.normal(size=7) + 1j*rng.normal(size=7))
        y2 = y0 + scale*(rng.normal(size=7) + 1j*rng.normal(size=7))
        new = []
        for S in sols:
            X = track(S, y0, y1, E)
            if X is None: continue
            X = track(X, y1, y2, E)
            if X is None: continue
            X = track(X, y2, y0, E)
            if X is None: continue
            if not any(np.linalg.norm(X - T) < 1e-6*(1+np.linalg.norm(T)) for T in sols + new):
                new.append(X)
        sols += new
        if verbose: print('loop', L, 'num sols', len(sols), flush=True)
        stale = 0 if new else stale + 1
        if stale >= 8: break
    return sols
if __name__ == '__main__':
    idx, k = sys.argv[1], sys.argv[2]
    V = np.load('m12_type%s_%s.npy' % (idx, k)); E = eval(open('m12_type%s_%s_E.txt' % (idx, k)).read())
    rng = np.random.default_rng(int(sys.argv[3]) if len(sys.argv) > 3 else 0)
    X, y = member_to_X(V, E, t0=0.3+0.7j)
    print('residual', np.linalg.norm(mob.F(X, y, E)), 'cond JX', np.linalg.cond(mob.J(X, y, E)))
    print('contacts', mob.contacts(X, y, E)[1])
    sols = monodromy(X, y, E, rng)
    flags = [mob.contacts(S, y, E)[1] for S in sols]
    print('type', idx, 'deg(ev) >=', len(sols), 'flags', flags)
