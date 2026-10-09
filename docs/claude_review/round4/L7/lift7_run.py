import numpy as np, sys
from lift7 import Sys
v0 = np.load('gen6_seed.npy')
seed = int(sys.argv[1]); iters = int(sys.argv[2])
rng = np.random.default_rng(seed); S = Sys(v0, rng)
u = np.concatenate([v0[S.free] + 0.05*(rng.normal(size=25)+1j*rng.normal(size=25)), rng.normal(size=2)+1j*rng.normal(size=2), rng.normal(size=18)+1j*rng.normal(size=18)])
n = len(u)
for it in range(iters):
    F = S.residual(u); nF = np.linalg.norm(F)
    if nF < 1e-12: break
    J = np.zeros((len(F), n), complex); h = 1e-7
    for k in range(n):
        du = np.zeros(n, complex); du[k] = h
        J[:, k] = (S.residual(u+du) - S.residual(u-du))/(2*h)
    du = np.linalg.lstsq(J, -F, rcond=None)[0]
    st = 1.0
    while st > 1e-12:
        F1 = S.residual(u + st*du)
        if np.all(np.isfinite(F1)) and np.linalg.norm(F1) < nF*(1-0.05*st): break
        st /= 2
    if st <= 1e-12: break
    u = u + st*du
sv = np.linalg.svd(J, compute_uv=False)
print('seed', seed, 'iters', it, 'resid', nF, 'J smallest sv', sv[-3:], flush=True)
np.save('lift7_run_%d.npy' % seed, u); np.save('lift7_runL_%d.npy' % seed, S.L)
