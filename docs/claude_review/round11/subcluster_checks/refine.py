# Gauss-Newton refinement of tangent-cone zeros and a robust (exhaustive, small-coefficient) degeneracy test.
import numpy as np, itertools
def resid_vec(systems, u, m_of):
    """stacked residual R(u) = [(I-QQ^T)(A_J u_J)^m / ||(A_J u_J)^m||]_J and Jacobian (w.r.t. u in R^M)."""
    Rs = []; Js = []
    M = len(u)
    for S in systems:
        J = list(S['J'] if isinstance(S, dict) else S.J); A = S['A'] if isinstance(S, dict) else S.A
        Q = S['Q'] if isinstance(S, dict) else S.Q; m = m_of(S)
        a = A @ u[J]; f = a**m; nf = np.linalg.norm(f)
        g = f - Q @ (Q.T @ f)
        df = (m * a**(m-1))[:, None] * A            # |A| x n
        dg = df - Q @ (Q.T @ df)
        # R = g/nf ; dR = dg/nf - g (f.df)/nf^3
        dR = dg/nf - np.outer(g, f @ df)/nf**3
        Jfull = np.zeros((len(a), M)); Jfull[:, J] = dR
        Rs.append(g/nf); Js.append(Jfull)
    return np.concatenate(Rs), np.vstack(Js)
def refine(systems, u, m_of, iters=60):
    M = len(u)
    u = u - u.mean(); u /= np.linalg.norm(u)
    for it in range(iters):
        R, Jm = resid_vec(systems, u, m_of)
        # restrict to tangent space of the sphere in R^M_0
        P = np.eye(M) - 1.0/M - np.outer(u, u)
        Jt = Jm @ P
        step, *_ = np.linalg.lstsq(Jt, -R, rcond=1e-12)
        u = u + P @ step; u -= u.mean(); u /= np.linalg.norm(u)
    R, _ = resid_vec(systems, u, m_of)
    return u, np.linalg.norm(R)
def min_relation(u, vmax):
    """min over 0 != v in Z^M_0 with |v|_inf <= vmax of |<v,u>|/|v|_1, and the minimiser (exhaustive)."""
    M = len(u); d = u[:-1] - u[-1]
    best = (np.inf, None)
    for v in itertools.product(range(-vmax, vmax+1), repeat=M-1):
        if not any(v): continue
        val = abs(np.dot(v, d)); l1 = sum(map(abs, v)) + abs(sum(v))
        if val/l1 < best[0]: best = (val/l1, list(v) + [-sum(v)])
    return best
