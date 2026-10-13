# Forward formulation: events P_e = (s_e, z_e); label a has (1,1)-curve c_a . (s z, s, z, 1) = 0
# through P_e for e containing a. Find members by Gauss-Newton; check nondegeneracy.
import numpy as np, itertools
def setup(E, rng):
    m = len(E); n = 7
    R = [rng.normal(size=4)+1j*rng.normal(size=4) for a in range(n)]  # normalization forms
    return R
def Fw(V, E, R):
    m = len(E); n = 7
    s = V[0:m]; z = V[m:2*m]; C = V[2*m:].reshape(n, 4)
    out = []
    for i, ev in enumerate(E):
        mon = np.array([s[i]*z[i], s[i], z[i], 1])
        for a in ev: out.append(C[a] @ mon)
    for a in range(n): out.append(R[a] @ C[a] - 1)
    return np.array(out)
def Jw(V, E, R):
    m = len(E); n = 7
    s = V[0:m]; z = V[m:2*m]; C = V[2*m:].reshape(n, 4)
    rows = []
    N = len(V)
    for i, ev in enumerate(E):
        mon = np.array([s[i]*z[i], s[i], z[i], 1]); ds = np.array([z[i], 1, 0, 0]); dz = np.array([s[i], 0, 1, 0])
        for a in ev:
            r = np.zeros(N, dtype=complex); r[i] = C[a] @ ds; r[m+i] = C[a] @ dz
            r[2*m+4*a:2*m+4*a+4] = mon; rows.append(r)
    for a in range(n):
        r = np.zeros(N, dtype=complex); r[2*m+4*a:2*m+4*a+4] = R[a]; rows.append(r)
    return np.array(rows)
def gn(V, E, R, maxit=100, tol=1e-13):
    for it in range(maxit):
        f = Fw(V, E, R)
        if np.linalg.norm(f) < tol: return V, True
        dV = np.linalg.lstsq(Jw(V, E, R), -f, rcond=None)[0]
        V = V + dV
        if not np.all(np.isfinite(V)) or np.linalg.norm(V) > 1e8: return V, False
    return V, np.linalg.norm(Fw(V, E, R)) < 1e-10
def to_mobius(c):
    # c0 s z + c1 s + c2 z + c3 = 0  ->  z = -(c1 s + c3)/(c0 s + c2): matrix [[-c1,-c3],[c0,c2]]
    return np.array([[-c[1], -c[3]], [c[0], c[2]]])
def check(V, E, tol=1e-6):
    m = len(E); n = 7
    s = V[0:m]; z = V[m:2*m]; C = V[2*m:].reshape(n, 4)
    g = [to_mobius(C[a]) for a in range(n)]
    for a in range(n):
        if abs(np.linalg.det(g[a])) < tol*np.linalg.norm(g[a])**2: return None, 'singular'
    times = [(s[i], frozenset(ev)) for i, ev in enumerate(E)]
    for a in range(n):
        for b in range(a+1, n):
            # fixed points of g_b^{-1} g_a : g_a(t) = g_b(t)
            A, B = g[a], g[b]
            # (A00 t + A01)(B10 t + B11) - (B00 t + B01)(A10 t + A11) = 0
            c2 = A[0,0]*B[1,0] - B[0,0]*A[1,0]; c1 = A[0,0]*B[1,1] + A[0,1]*B[1,0] - B[0,0]*A[1,1] - B[0,1]*A[1,0]
            c0 = A[0,1]*B[1,1] - B[0,1]*A[1,1]
            if abs(c2) < 1e-10 and abs(c1) < 1e-10 and abs(c0) < 1e-10: return None, 'equal'
            roots = np.roots([c2, c1, c0]) if abs(c2) > 1e-12 else np.roots([c1, c0])
            for t in roots:
                if not any(abs(t - u) < tol*(1+abs(u)) and {a, b} <= S for (u, S) in times[:m]):
                    times.append((t, frozenset((a, b))))
    us = [t for t, S in times]
    for i in range(len(us)):
        for j in range(i+1, len(us)):
            if abs(us[i]-us[j]) < tol*(1+abs(us[i])): return None, 'coincide'
    for (t, S) in times:
        pos = [(g[a][0,0]*t + g[a][0,1])/(g[a][1,0]*t + g[a][1,1]) for a in range(n)]
        cl = [pos[a] for a in S]
        if max(abs(c - cl[0]) for c in cl) > 1e-5*(1+abs(cl[0])): return None, 'notcluster'
        others = [pos[a] for a in range(n) if a not in S] + [cl[0]]
        for i in range(len(others)):
            for j in range(i+1, len(others)):
                if abs(others[i]-others[j]) < tol*(1+abs(others[i])): return None, 'extra'
    return (g, times), 'ok'
