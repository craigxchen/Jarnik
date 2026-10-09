"""NUMERICAL control: same LM solver on the k=6 C5 ansatz (cubics), where nondegenerate solutions exist."""
import numpy as np, itertools, sys
L = [1, 2, 3, 4, 5]
EDGES = [(1, 2), (2, 3), (3, 4), (4, 5), (1, 5)]
TRIP = [(1, 2, 4), (2, 3, 5), (1, 3, 4), (2, 4, 5), (1, 3, 5)]
def nbr(i):
    return [e for e in EDGES if i in e]
def unpack(z):
    r_e = dict(zip(EDGES, z[0:5])); r = dict(zip(L, z[5:10])); c = dict(zip(L, np.concatenate([[1.0], z[10:14]])))
    sig = dict(zip(TRIP, z[14:19])); return r_e, r, c, sig
def u(i, s, r_e, r, c):
    v = c[i] * (s - r[i])
    for e in nbr(i): v = v * (s - r_e[e])
    return v
def F(z):
    r_e, r, c, sig = unpack(z); out = []
    for T in TRIP:
        a, b, cc = T; s = sig[T]
        ua, ub, uc = u(a, s, r_e, r, c), u(b, s, r_e, r, c), u(cc, s, r_e, r, c)
        out += [ua - ub, ua - uc]
    out += [r_e[(1, 2)], r_e[(2, 3)] - 1]
    return np.array(out)
def J(z, h=1e-7):
    f0 = F(z); Jm = np.zeros((len(f0), len(z)), dtype=complex)
    for k in range(len(z)):
        dz = np.zeros(len(z), dtype=complex); dz[k] = h; Jm[:, k] = (F(z + dz) - f0) / h
    return Jm
def degeneracy(z):
    r_e, r, c, sig = unpack(z)
    vals = list(r_e.values()) + list(r.values()) + list(sig.values())
    md = min(abs(x - y) for x, y in itertools.combinations(vals, 2))
    cs = list(c.values())
    mc = min(min(abs(x - y) for x, y in itertools.combinations(cs, 2)), min(abs(x) for x in cs))
    return md, mc
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
good = 0; ntr = 30
for trial in range(ntr):
    z = (rng.normal(size=19) + 1j * rng.normal(size=19)) * 2; lam = 1e-3
    for it in range(400):
        f = F(z); nf = np.linalg.norm(f)
        if nf < 1e-12: break
        Jm = J(z); A = Jm.conj().T @ Jm + lam * np.eye(len(z))
        z2 = z + np.linalg.solve(A, -Jm.conj().T @ f)
        if np.linalg.norm(F(z2)) < nf: z = z2; lam = max(lam / 3, 1e-12)
        else:
            lam *= 4
            if lam > 1e8: break
    md, mc = degeneracy(z); res = np.linalg.norm(F(z))
    ok = res < 1e-9 and md > 1e-3 and mc > 1e-3
    good += ok
print('k=6 control: nondegenerate solutions found in %d of %d trials' % (good, ntr))
