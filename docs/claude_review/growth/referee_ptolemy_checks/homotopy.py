r"""NUMERICAL exploration (not a proof): balanced degree-one rational curves in Mbar_{0,m+1}
by Newton-homotopy path tracking in the X-chart.

Unknowns z: coefficients of X_2..X_m below the (fixed, random) leading coefficients, and the
roots rho_T (3<=|T|<=m-1) except a fixed random subset making the system square.
F(z) = (X_a(rho_T) - X_b(rho_T))_{T, consecutive a<b in T}.
Start: random z0 (all roots distinct).  Homotopy H(z,t) = F(z) - phi(t) F(z0), phi(t) = (1-t)(1+(gamma-1)t), t: 0 -> 1,
gamma a random unit complex number (complex path, gamma trick).  Predictor-corrector with adaptive steps.
Endpoints are checked for nondegeneracy (all 2^m-m-2 boundary roots distinct, leads distinct).
Usage: python3 homotopy.py m trials [seed]
"""
import sys, itertools
import numpy as np
from balanced_newton import setup, X, dX, pair_roots, verify

class Sys:
    def __init__(self, m, rng):
        self.m = m
        rows, big, D, conds = setup(m)
        self.rows, self.big, self.D, self.conds = rows, big, D, conds
        nb = len(big); E = len(conds); ncoef = (m - 1) * D
        nfix = ncoef + nb - E
        assert 0 <= nfix <= nb
        c = lambda n: (rng.normal(size=n) + 1j * rng.normal(size=n)) / np.sqrt(2)
        self.leads = dict(zip(rows[1:], c(m - 1)))
        order = list(range(nb)); rng.shuffle(order)
        fixed = set(order[:nfix])
        self.rho0 = {T: v for T, v in zip(big, c(nb) * 2)}
        self.free = [T for n, T in enumerate(big) if n not in fixed]
        self.fixedvals = {T: self.rho0[T] for n, T in enumerate(big) if n in fixed}
        self.ncoef = ncoef; self.E = E
        self.z0 = np.concatenate([c(ncoef), np.array([self.rho0[T] for T in self.free])])
        self.rpos = {T: ncoef + n for n, T in enumerate(self.free)}
    def unpack(self, z):
        D = self.D
        coef = {i: np.concatenate([[self.leads[i]], z[(i - 2) * D:(i - 1) * D]]) for i in self.rows[1:]}
        rho = dict(self.fixedvals)
        for T in self.free:
            rho[T] = z[self.rpos[T]]
        return coef, rho
    def F(self, z):
        coef, rho = self.unpack(z)
        return np.array([X(coef, a, rho[T]) - X(coef, b, rho[T]) for T, a, b in self.conds])
    def J(self, z):
        coef, rho = self.unpack(z)
        D = self.D
        Jm = np.zeros((self.E, len(z)), dtype=complex)
        for e, (T, a, b) in enumerate(self.conds):
            x = rho[T]
            pw = x ** np.arange(D - 1, -1, -1)
            for i, s in ((a, 1), (b, -1)):
                if i > 1:
                    Jm[e, (i - 2) * D:(i - 1) * D] += s * pw
            if T in self.rpos:
                Jm[e, self.rpos[T]] = dX(coef, a, x) - dX(coef, b, x)
        return Jm

def track(S, rng, maxsteps=4000):
    gamma = np.exp(2j * np.pi * rng.random())
    F0 = S.F(S.z0)
    z = S.z0.copy(); t = 0.0; h = 0.01
    phi = lambda t: (1 - t) * (1 + (gamma - 1) * t)
    dphi = lambda t: -(1 + (gamma - 1) * t) + (1 - t) * (gamma - 1)
    H = lambda z, t: S.F(z) - phi(t) * F0
    steps = 0
    while t < 1.0 and steps < maxsteps:
        steps += 1
        h = min(h, 1.0 - t)
        Jm = S.J(z)
        try:
            dzdt = np.linalg.solve(Jm, dphi(t) * F0)     # d/dt H = 0  ->  J dz/dt - phi'(t) F0 = 0
        except np.linalg.LinAlgError:
            return None, 'singular'
        zp = z + h * dzdt
        tn = t + h
        ok = False
        for it in range(6):
            r = H(zp, tn)
            try:
                dz = np.linalg.solve(S.J(zp), -r)
            except np.linalg.LinAlgError:
                break
            zp = zp + dz
            if np.linalg.norm(dz) < 1e-9 * (1 + np.linalg.norm(zp)):
                ok = True
                break
        if ok and np.linalg.norm(H(zp, tn)) < 1e-8 * (1 + np.linalg.norm(F0)):
            z, t = zp, tn
            h = min(h * 1.5, 0.05)
        else:
            h /= 2
            if h < 1e-10:
                return None, f'step underflow at t={t:.6f}'
    # polish
    for it in range(30):
        f = S.F(z)
        if np.linalg.norm(f) < 1e-13:
            break
        try:
            z = z + np.linalg.solve(S.J(z), -f)
        except np.linalg.LinAlgError:
            break
    return z, f't=1 reached in {steps} steps'

if __name__ == "__main__":
    m = int(sys.argv[1]); trials = int(sys.argv[2]); seed = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    rng = np.random.default_rng(seed)
    good = 0
    for tr in range(trials):
        S = Sys(m, rng)
        z, msg = track(S, rng)
        if z is None:
            print(f"trial {tr}: {msg}", flush=True); continue
        nr = np.linalg.norm(S.F(z))
        coef, rho = S.unpack(z)
        pr = pair_roots(m, S.rows, S.big, coef, rho)
        allrho = dict(rho); allrho.update(pr)
        vals = list(allrho.values())
        sep = min(abs(a - b) for a, b in itertools.combinations(vals, 2))
        leads = [0.0] + [coef[i][0] for i in S.rows[1:]]
        lsep = min(abs(x - y) for x, y in itertools.combinations(leads, 2))
        ok = nr < 1e-9 and sep > 1e-5 and lsep > 1e-5
        w = verify(m, S.rows, coef, allrho) if ok else float('nan')
        good += ok
        print(f"trial {tr}: {msg}; residual {nr:.1e}, min sep {sep:.2e}, identity err {w:.1e} -> {'NONDEGENERATE' if ok else 'degenerate'}", flush=True)
        if ok:
            np.save(f"sol_m{m}_seed{seed}_tr{tr}.npy", np.array([z], dtype=object), allow_pickle=True)
    print("nondegenerate:", good, "/", trials)
