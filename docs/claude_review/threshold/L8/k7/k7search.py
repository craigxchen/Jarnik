"""NUMERICAL (evidence only). k=7 balanced rational curves in the infinity-unmarked chart:
u_0 = 0, u_i(s) = c_i (s - r_i) prod_{j != i} (s - r_ij)  (i=1..6, r_ij = r_ji: 0,i,j collide),
and for every triple T={a,b,c} of {1..6}: u_a(sig_T) = u_b(sig_T) = u_c(sig_T)."""
import numpy as np, itertools, sys
L = range(1, 7)
EDGES = list(itertools.combinations(L, 2))
TRIP = list(itertools.combinations(L, 3))
def unpack(z):
    r_e = dict(zip(EDGES, z[0:15])); r = dict(zip(L, z[15:21])); c = dict(zip(L, np.concatenate([[1.0], z[21:26]])))
    sig = dict(zip(TRIP, z[26:46]))
    return r_e, r, c, sig
def u(i, s, r_e, r, c):
    v = c[i] * (s - r[i])
    for j in L:
        if j != i:
            v = v * (s - r_e[tuple(sorted((i, j)))])
    return v
def F(z):
    r_e, r, c, sig = unpack(z)
    out = []
    for T in TRIP:
        a, b, cc = T; s = sig[T]
        ua, ub, uc = u(a, s, r_e, r, c), u(b, s, r_e, r, c), u(cc, s, r_e, r, c)
        out += [ua - ub, ua - uc]
    # affine normalization: r_12 = 0, r_13 = 1
    out += [r_e[(1, 2)], r_e[(1, 3)] - 1]
    return np.array(out)
def J(z, h=1e-7):
    f0 = F(z); Jm = np.zeros((len(f0), len(z)), dtype=complex)
    for k in range(len(z)):
        dz = np.zeros(len(z), dtype=complex); dz[k] = h
        Jm[:, k] = (F(z + dz) - f0) / h
    return Jm
def degeneracy(z):
    r_e, r, c, sig = unpack(z)
    vals = list(r_e.values()) + list(r.values()) + list(sig.values())
    md = min(abs(x - y) for x, y in itertools.combinations(vals, 2))
    cs = list(c.values())
    mc = min(min(abs(x - y) for x, y in itertools.combinations(cs, 2)), min(abs(x) for x in cs))
    return md, mc, max(abs(x) for x in vals)
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
ntr = int(sys.argv[2]) if len(sys.argv) > 2 else 40
good = 0
for trial in range(ntr):
    z = (rng.normal(size=46) + 1j * rng.normal(size=46)) * 2
    lam = 1e-3
    for it in range(400):
        f = F(z); nf = np.linalg.norm(f)
        if nf < 1e-12: break
        Jm = J(z)
        # Levenberg-Marquardt minimum-norm step
        A = Jm.conj().T @ Jm + lam * np.eye(len(z))
        step = np.linalg.solve(A, -Jm.conj().T @ f)
        z2 = z + step
        if np.linalg.norm(F(z2)) < nf:
            z = z2; lam = max(lam / 3, 1e-12)
        else:
            lam *= 4
            if lam > 1e8: break
    md, mc, mx = degeneracy(z)
    res = np.linalg.norm(F(z))
    tag = 'NONDEGENERATE?' if res < 1e-9 and md > 1e-3 and mc > 1e-3 and mx < 1e4 else ''
    if res < 1e-9:
        good += tag != ''
    print('trial %3d  |F|=%.1e  min sep=%.1e  min c-sep=%.1e  max|val|=%.1e  %s' % (trial, res, md, mc, mx, tag), flush=True)
print('nondegenerate candidates:', good)
