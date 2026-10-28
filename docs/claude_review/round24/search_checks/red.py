# reduced system: unknowns r_T, |T|>=2 (11 complex), gauge Re r_1234 = 0
import itertools, numpy as np
BIG = [frozenset(c) for k in (2,3,4) for c in itertools.combinations(range(1,5),k)]  # 11; last is 1234
BIDX = {T:i for i,T in enumerate(BIG)}
SMEM = [[BIDX[T] for T in BIG if j in T] for j in range(1,5)]
W = np.poly1d([1,0,-1])**4
Wc = W.coeffs.astype(float)  # 9 coeffs high->low

def unpack(x):
    r = np.zeros(11, dtype=complex)
    r[:10] = x[:10] + 1j*x[10:20]
    r[10] = 1j*x[20]
    return r

def blocks(r, j):
    P = np.array([1.0+0j])
    for k in SMEM[j]:
        P = np.convolve(P, [1.0, -r[k]])
    # 9-coefficient vectors: Pi_j (deg 7, padded), t*Pi_j (deg 8)
    P9 = np.concatenate([[0], P]); tP = np.concatenate([P, [0]])
    return P9, tP

def solve_j(r, j):
    P9, tP = blocks(r, j)
    V = np.stack([P9.real, P9.imag, tP.real, tP.imag], axis=1)  # 9x4
    coef, *_ = np.linalg.lstsq(V, Wc, rcond=None)
    res = Wc - V@coef
    # Im((lam t + mu) Pi) = Im(lam) Re tP + Re(lam) Im tP + Im(mu) Re P + Re(mu) Im P = W
    mu = coef[1] + 1j*coef[0]; lam = coef[3] + 1j*coef[2]
    return res, lam, mu

def residual(x):
    r = unpack(x)
    return np.concatenate([solve_j(r,j)[0] for j in range(4)])

def full_roots(x):
    """return dict T -> r_T for all 15 T and the g_j"""
    r = unpack(x); out = {T: r[BIDX[T]] for T in BIG}; gs = []
    for j in range(4):
        res, lam, mu = solve_j(r, j)
        out[frozenset([j+1])] = -mu/lam
        P9, tP = blocks(r, j)
        gs.append(lam*tP + mu*P9)
    return out, gs

def diagnostics(x):
    R, gs = full_roots(x)
    vals = np.array(list(R.values()))
    imin = np.min(np.abs(vals.imag))
    sep = min(abs(a-b) for a,b in itertools.combinations(vals,2))
    csep = min(abs(a-np.conj(b)) for a in vals for b in vals)
    tt = 0.3+0.0j
    us = [np.polyval(g,tt)/np.conj(np.polyval(g,tt)) for g in gs]
    du = min([abs(a-b) for a,b in itertools.combinations(us,2)]+[abs(u-1) for u in us])
    # check g_j = a_j + i W with a_j real: Im g_j coefficients == W
    chk = max(np.max(np.abs(g.imag - Wc)) for g in gs)
    return dict(imin=imin, sep=sep, csep=csep, du=du, rmax=np.max(np.abs(vals)), chk=chk)
