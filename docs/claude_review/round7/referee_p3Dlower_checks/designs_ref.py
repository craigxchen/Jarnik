"""Referee's independent recomputation of the p3Dlower §6 evidence table, with CORRECT per-unit bookkeeping
(unit s defined by <v,ell> = s*iota mod m_q, iota = log i at that q), and, for comparison, the 'raw class'
bookkeeping u = <v,ell>/(m_q/4) used by the author's designs.py (which mixes units i and -i across q).
Designs as in the write-up: RAND (uniform in parity class), COH (zeta_j = smallest split primes > X),
TWIST (Def. 4.1), GEN (genuine pi_j, p_j ~ M^3 drawn with random.seed(1) exactly as the author's script).
usage: python3 designs_ref.py M kmax"""
import sys, math, random
import numpy as np
from sympy import primerange, isprime, legendre_symbol
from tq import T, two_sq, D
M = int(sys.argv[1]); KMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 2
r = M; X = 3 * M; QMIN = 0.5
random.seed(1)
lo = int(M ** 3.0); P = set()
while len(P) < r:
    c = random.randrange(lo, 2 * lo)
    if c % 4 == 1 and c > X and isprime(c): P.add(c)
P = sorted(P); pis = [two_sq(p) for p in P]
split = [p for p in primerange(5, 10**6) if p % 4 == 1]
zs = [two_sq(p) for p in [p for p in split if p > X][:r]]
rng = np.random.default_rng(12345)
qs = list(primerange(3, X + 1))
labs = {k: [] for k in ('RAND', 'COH', 'TWIST', 'GEN')}
ms, iotas = [], []
for q in qs:
    t = T(q); m = t.m; ms.append(m); iotas.append(t.iota)
    nu = np.array([1 if legendre_symbol(p % q, q) == -1 else 0 for p in P])
    labs['RAND'].append((2 * rng.integers(0, m // 2, size=r) + nu) % m)
    coh = np.array([t.L(a, b) for a, b in zs]); labs['COH'].append(coh)
    delta = (coh % 2) ^ nu
    if q % 8 in (3, 5): tw = t.iota
    else:
        ell = next(l for l in split if legendre_symbol(l % q, q) == -1); tw = t.L(*two_sq(ell))
    labs['TWIST'].append((coh + delta * tw) % m)
    gen = np.array([t.L(a, b) for a, b in pis]); labs['GEN'].append(gen)
    assert np.all(labs['TWIST'][-1] % 2 == nu) and np.all(gen % 2 == nu) and np.all(labs['RAND'][-1] % 2 == nu)
ms = np.array(ms); io = np.array(iotas); qa = np.array(qs)
logq = np.log(qa.astype(float)) * (qa >= M ** QMIN); lM = math.log(M)
sfac = np.where(io == ms // 4, 1, 3)        # u -> s = u*sfac mod 4
def demand(V, correct=True):
    sh = (-1,) + (1,) * (V.ndim - 1)
    Vm = V % ms.reshape(sh); qq = (ms // 4).reshape(sh)
    hit = (Vm % qq == 0); u = Vm // qq
    s = (u * sfac.reshape(sh)) % 4 if correct else u
    lq = logq.reshape(sh)
    return np.max(np.stack([((hit & (s == k)) * lq).sum(0) for k in range(4)]), axis=0)
print(f"M={M} r={r} X={X} #moduli={len(qs)}; budget w(v)/2 ~ 1.5 k log M (+0.35 per k)")
for name, Ls in labs.items():
    Lm = np.stack(Ls); out = []
    for corr in (True, False):
        res = [demand(Lm, corr).max() / lM]
        b2 = demand(2 * Lm, corr).max()
        for sg in (1, -1):
            d = demand(Lm[:, :, None] + sg * Lm[:, None, :], corr); np.fill_diagonal(d, 0); b2 = max(b2, d.max())
        res.append(b2 / lM)
        if KMAX >= 3:
            b3 = demand(3 * Lm, corr).max()
            for sg in (1, -1):
                d = demand(2 * Lm[:, :, None] + sg * Lm[:, None, :], corr); np.fill_diagonal(d, 0); b3 = max(b3, d.max())
            for j in range(r):
                for s1 in (1, -1):
                    for s2 in (1, -1):
                        d = demand(Lm[:, j][:, None, None] + s1 * Lm[:, :, None] + s2 * Lm[:, None, :], corr)
                        d[j, :] = 0; d[:, j] = 0; np.fill_diagonal(d, 0); b3 = max(b3, d.max())
            res.append(b3 / lM)
        out.append(res)
    f = lambda xs: ", ".join(f"{x:5.2f}" for x in xs)
    print(f"  {name:5s} correct per-unit: k=1..{KMAX}: {f(out[0])}    raw-class (author's bookkeeping): {f(out[1])}")
