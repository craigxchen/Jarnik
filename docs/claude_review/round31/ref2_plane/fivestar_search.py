"""Vectorised Gauss-Newton search (min-norm steps) for real/complex y on the 10 edges of K5 with
constant star sums and constant reciprocal star sums, normalised sum y^2 = 10.  Real: expect only
degenerate solutions (some |y_e| equal).  Complex: expect nondegenerate solutions."""
import numpy as np, sys
edges = [(i, j) for i in range(5) for j in range(i+1, 5)]
N = np.zeros((5, 10))
for k, (i, j) in enumerate(edges):
    N[i, k] = N[j, k] = 1
Dm = N[:4] - N[4]          # 4x10 differences of star rows

def run(Y, iters=80):
    for it in range(iters):
        F = np.concatenate([Y @ Dm.T, (1/Y) @ Dm.T, (np.sum(Y*Y, axis=1) - 10)[:, None]], axis=1)  # B x 9
        J1 = np.broadcast_to(Dm, (Y.shape[0], 4, 10))
        J2 = Dm[None, :, :] * (-1/Y**2)[:, None, :]
        J3 = (2*Y)[:, None, :]
        J = np.concatenate([J1, J2, J3], axis=1)    # B x 9 x 10
        step = np.einsum('bij,bj->bi', np.linalg.pinv(J), F)
        Y = Y - step
    F = np.concatenate([Y @ Dm.T, (1/Y) @ Dm.T, (np.sum(Y*Y, axis=1) - 10)[:, None]], axis=1)
    return Y, np.linalg.norm(F, axis=1)

rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 7)
B = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
with np.errstate(all='ignore'):
    Y0 = rng.normal(size=(B, 10)) * rng.choice([0.3, 1.0, 3.0], size=(B, 1)) + rng.normal(size=(B,1))
    Y, res = run(Y0)
ok = np.isfinite(res) & (res < 1e-9) & (np.min(np.abs(Y), axis=1) > 1e-6)
Yr = Y[ok]
a = np.sort(np.abs(Yr), axis=1)
gap = np.min(np.diff(a, axis=1), axis=1) / np.max(a, axis=1)
print("REAL: starts", B, "converged", ok.sum())
if ok.sum():
    print("REAL: max over solutions of (min gap between |y_e|)/max|y|: %.3e" % gap.max())
    print("REAL: quantiles of gap:", np.quantile(gap, [0.5, 0.99, 1.0]))
    print("REAL: example solution:", np.round(Yr[np.argmax(gap)], 5))
sys.stdout.flush()

# complex
Bc = 4000
with np.errstate(all='ignore'):
    Yc = (rng.normal(size=(Bc, 10)) + 1j*rng.normal(size=(Bc, 10)))
    for it in range(120):
        F = np.concatenate([Yc @ Dm.T, (1/Yc) @ Dm.T, (np.sum(Yc*Yc, axis=1) - 10)[:, None]], axis=1)
        J = np.concatenate([np.broadcast_to(Dm, (Bc, 4, 10)).astype(complex),
                            Dm[None] * (-1/Yc**2)[:, None, :], (2*Yc)[:, None, :]], axis=1)
        Yc = Yc - np.einsum('bij,bj->bi', np.linalg.pinv(J), F)
    F = np.concatenate([Yc @ Dm.T, (1/Yc) @ Dm.T, (np.sum(Yc*Yc, axis=1) - 10)[:, None]], axis=1)
    rc = np.linalg.norm(F, axis=1)
okc = np.isfinite(rc) & (rc < 1e-9)
Z = Yc[okc]
y2 = Z**2
sep = np.array([min(abs(r[i]-r[j]) for i in range(10) for j in range(i+1, 10)) for r in y2]) if len(Z) else np.array([])
nd = (sep > 1e-3) & (np.min(np.abs(Z), axis=1) > 1e-3) if len(Z) else np.array([], bool)
print("COMPLEX: converged", okc.sum(), "nondegenerate", nd.sum())
if nd.sum():
    Zn = Z[nd]
    # how close to real can a nondegenerate complex solution be, after the best common complex rotation?
    # (solutions are only defined up to scaling y -> c y with c^2 real-normalised; rotate by +-1, +-i)
    best = []
    for z in Zn:
        m = min(np.max(np.abs((w*z).imag)) / np.max(np.abs(z)) for w in [1, 1j])
        best.append(m)
    print("COMPLEX: min over nondegenerate solutions of max|Im| / max|y| (up to y->iy): %.3e" % min(best))
    # check the lemma's prediction: s and r (star sums, reciprocal sums)
    z = Zn[0]
    print("COMPLEX example star sums", np.round(N @ z, 6), "recip", np.round(N @ (1/z), 6))
