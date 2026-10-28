# twisted two-branch ansatz: u_j = G_j/conj(G_j), G_j = e^{i th_j} prod_{T contains j} (t - r_T), T subset {1..4}
import itertools, numpy as np
SUBS = [frozenset(c) for k in range(1,5) for c in itertools.combinations(range(1,5),k)]  # 15
IDX = {T:i for i,T in enumerate(SUBS)}
MEM = [[IDX[T] for T in SUBS if j in T] for j in range(1,5)]

def taylor_at(roots, phase, t0, K=4):
    """first K Taylor coeffs at t0 of phase*prod(t - r)"""
    c = np.zeros(K, dtype=complex); c[0] = phase
    for r in roots:
        a = t0 - r   # (t - r) = a + h
        new = c*a
        new[1:] += c[:-1]
        c = new
    return c

def residual(x, t0s=(1.0,-1.0)):
    a = x[:15]; b = x[15:30]; th = x[30:34]
    r = a + 1j*b
    res = []
    for j in range(4):
        ph = np.exp(1j*th[j])
        for t0 in t0s:
            c = taylor_at(r[MEM[j]], ph, t0)
            res.extend(c.imag)
    return np.array(res)

def Gpoly(x, j):
    r = x[:15] + 1j*x[15:30]
    p = np.array([np.exp(1j*x[30+j])])
    for k in MEM[j]:
        p = np.convolve(p, np.array([1.0, -r[k]]))
    return p
