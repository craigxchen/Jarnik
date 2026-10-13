# Mobius curves through a fixed point x = (y_0..y_6): inverse maps G_a(z) = (p_a z + q_a)/(z - y_a),
# gauge G_0 = 1/(z - y_0). Events e = (a,b,c): z_e with G_a(z_e)=G_b(z_e)=G_c(z_e).
# Unknown vector X = (p_1..p_6, q_1..q_6, z_0..z_{m-1}).
import numpy as np
rng = np.random.default_rng()

def unpack(X, n=7):
    p = np.concatenate([[0], X[0:n-1]]); q = np.concatenate([[1], X[n-1:2*n-2]]); z = X[2*n-2:]
    return p, q, z

def eqlist(E):
    L = []
    for i, ev in enumerate(E):
        if len(ev) == 3:
            a, b, c = ev; L += [(i, a, b), (i, b, c)]
        else:  # 4-event
            a, b, c, d = ev; L += [(i, a, b), (i, b, c), (i, c, d)]
    return L

def F(X, y, E):
    p, q, z = unpack(X)
    out = []
    for (i, a, b) in eqlist(E):
        zz = z[i]
        out.append((p[a]*zz + q[a])*(zz - y[b]) - (p[b]*zz + q[b])*(zz - y[a]))
    return np.array(out)

def J(X, y, E, wrt_y=False):
    n = 7; p, q, z = unpack(X); m = len(E)
    L = eqlist(E)
    JX = np.zeros((len(L), len(X)), dtype=complex); Jy = np.zeros((len(L), n), dtype=complex)
    for r, (i, a, b) in enumerate(L):
        zz = z[i]
        if a > 0:
            JX[r, a-1] += zz*(zz - y[b]); JX[r, n-1+a-1] += (zz - y[b])
        if b > 0:
            JX[r, b-1] -= zz*(zz - y[a]); JX[r, n-1+b-1] -= (zz - y[a])
        JX[r, 2*n-2+i] = p[a]*(zz - y[b]) + (p[a]*zz + q[a]) - p[b]*(zz - y[a]) - (p[b]*zz + q[b])
        Jy[r, b] += -(p[a]*zz + q[a]); Jy[r, a] += (p[b]*zz + q[b])
    return (JX, Jy) if wrt_y else JX

def newton(X, y, E, tol=1e-12, maxit=60):
    for it in range(maxit):
        f = F(X, y, E)
        if not np.all(np.isfinite(f)): return X, False
        if np.linalg.norm(f) < tol*(1+np.linalg.norm(X)): return X, True
        Jm = J(X, y, E)
        try: dX = np.linalg.lstsq(Jm, -f, rcond=None)[0]
        except Exception: return X, False
        X = X + dX
        if np.linalg.norm(X) > 1e8: return X, False
    f = F(X, y, E)
    return X, np.linalg.norm(f) < 1e-9*(1+np.linalg.norm(X))

def Ginv(p, q, y, a, u):
    # forward position of label a at parameter u: solve (p z + q) = u (z - y)
    return (q[a] + u*y[a])/(u - p[a])

def contacts(X, y, E, tol=1e-6):
    """Return list of contact parameter values u with the set of labels colliding, and a
    nondegeneracy flag. Checks: G_a nonconstant, every event a genuine triple at a time u_e
    distinct from all other contact times, pure pair collisions distinct."""
    n = 7; p, q, z = unpack(X)
    for a in range(n):
        if abs(p[a]*y[a] + q[a]) < tol*(1+abs(p[a])+abs(q[a])): return None, 'const'
    times = []  # (u, set)
    for i, ev in enumerate(E):
        u = (p[ev[0]]*z[i] + q[ev[0]])/(z[i] - y[ev[0]])
        times.append((u, frozenset(ev)))
    # pair collisions: roots of F_ab
    for a in range(n):
        for b in range(a+1, n):
            c2 = p[a] - p[b]; c1 = q[a] - q[b] - p[a]*y[b] + p[b]*y[a]; c0 = q[b]*y[a] - q[a]*y[b]
            if abs(c2) < 1e-12: return None, 'deg'
            roots = np.roots([c2, c1, c0])
            for zz in roots:
                u = (p[a]*zz + q[a])/(zz - y[a])
                # is this an event time?
                if not any(abs(u - t) < tol*(1+abs(t)) and {a, b} <= s for (t, s) in times[:len(E)]):
                    times.append((u, frozenset((a, b))))
    # distinctness of all contact times, and genericity at each contact
    us = [t for t, s in times]
    for i in range(len(us)):
        for j in range(i+1, len(us)):
            if abs(us[i]-us[j]) < tol*(1+abs(us[i])): return None, 'coincide'
    # at each contact time the labels outside the set are distinct from each other and the cluster
    for (u, s) in times:
        pos = [Ginv(p, q, y, a, u) for a in range(n)]
        cl = [pos[a] for a in s]
        if max(abs(c - cl[0]) for c in cl) > 1e-5*(1+abs(cl[0])): return None, 'notcluster'
        others = [pos[a] for a in range(n) if a not in s] + [cl[0]]
        for i in range(len(others)):
            for j in range(i+1, len(others)):
                if abs(others[i]-others[j]) < tol*(1+abs(others[i])): return None, 'extra'
    return times, 'ok'

def classvec(times):
    pair = sum(1 for t, s in times if len(s) == 2 or len(s) == 5)
    trip = sum(1 for t, s in times if len(s) in (3, 4))
    return pair, trip
