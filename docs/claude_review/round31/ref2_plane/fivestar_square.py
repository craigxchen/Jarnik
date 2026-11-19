"""Square real system: 8 star equations + y_0 = 1 + y_1 = a, for a on a grid in (-3,3) minus {+-1, 0}.
Newton from many random real starts.  A real nondegenerate solution curve would show up for an interval of a."""
import numpy as np, sys
edges = [(i, j) for i in range(5) for j in range(i+1, 5)]
N = np.zeros((5, 10))
for k, (i, j) in enumerate(edges):
    N[i, k] = N[j, k] = 1
Dm = N[:4] - N[4]
rng = np.random.default_rng(11)
avals = np.concatenate([np.linspace(-3, 3, 61)])
avals = avals[(np.abs(np.abs(avals) - 1) > 1e-6) & (np.abs(avals) > 1e-6)]
B = 2000
tot_conv = 0; best = 0; bestY = None; perm_best = None
for perm_try in range(4):           # also vary which two edges are pinned
    e0, e1 = rng.choice(10, 2, replace=False)
    for a in avals:
        with np.errstate(all='ignore'):
            Y = rng.normal(size=(B, 10)) * 2
            Y[:, e0] = 1; Y[:, e1] = a
            for it in range(60):
                F = np.concatenate([Y @ Dm.T, (1/Y) @ Dm.T, (Y[:, e0]-1)[:, None], (Y[:, e1]-a)[:, None]], axis=1)
                E0 = np.zeros((B, 1, 10)); E0[:, 0, e0] = 1
                E1 = np.zeros((B, 1, 10)); E1[:, 0, e1] = 1
                J = np.concatenate([np.broadcast_to(Dm, (B, 4, 10)), Dm[None] * (-1/Y**2)[:, None, :], E0, E1], axis=1)
                JT = np.transpose(J, (0, 2, 1)); Y = Y - np.linalg.solve(JT @ J + 1e-12*np.eye(10), (JT @ F[..., None]))[..., 0]
            F = np.concatenate([Y @ Dm.T, (1/Y) @ Dm.T, (Y[:, e0]-1)[:, None], (Y[:, e1]-a)[:, None]], axis=1)
            r = np.linalg.norm(F, axis=1)
        ok = np.isfinite(r) & (r < 1e-9) & (np.min(np.abs(Y), axis=1) > 1e-6)
        tot_conv += ok.sum()
        if ok.sum():
            A = np.sort(np.abs(Y[ok]), axis=1)
            g = np.min(np.diff(A, axis=1), axis=1) / A[:, -1]
            if g.max() > best:
                best = g.max(); bestY = Y[ok][np.argmax(g)]; perm_best = (e0, e1, a)
    print("pin", (e0, e1), "cumulative converged", tot_conv, "best gap so far %.3e" % best, flush=True)
print("REAL square system: total converged", tot_conv, "max normalised min-gap %.3e" % best)
print("best solution", np.round(bestY, 6) if bestY is not None else None, perm_best)
