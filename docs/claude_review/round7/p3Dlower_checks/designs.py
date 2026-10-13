"""EVIDENCE ONLY. Monomial collision demands (level 1, all four units) for sparse v of norm k<=3, at odd primes
q <= X, for four label designs on r columns:
  RAND : uniform labels in the parity class nu_j(q) = [(p_j/q) = -1]           (in P3-all*; sharp.md random fibre type)
  COH  : ell_j = log(z_j/conj z_j), z_j small Gaussian primes                   (NOT in class: wrong parity)
  TWIST: COH + delta_j(q) * log(t_q), t_q a fixed non-square of T_q              (in class; Sec. 5 of p3Dlower.md)
  GEN  : ell_j = log(pi_j/conj pi_j), Norm pi_j = p_j                            (in class; genuine)
demand_s(v) = sum of log q over q with <v,ell(q)> = log(i^s) mod m_q (per unit s); report max_s.  Reports max over v of demand/log M, per norm k.
usage: python3 designs.py M r C seed"""
import sys, math, random
import numpy as np
from sympy import primerange, isprime, legendre_symbol
from tgroup import Tq, gauss_prime

M = int(sys.argv[1]); r = int(sys.argv[2]); Cexp = float(sys.argv[3]); seed = int(sys.argv[4])
KMAX = int(sys.argv[5]) if len(sys.argv) > 5 else 3
QMIN = float(sys.argv[6]) if len(sys.argv) > 6 else 0.0   # only moduli q >= M^QMIN count
random.seed(seed); np.random.seed(seed)
X = 3 * M
small = [p for p in primerange(5, 10**6) if p % 4 == 1]
zsm = [p for p in small if p > 3 * M][:r]
zeta = [gauss_prime(p) for p in zsm]                     # COH columns: smallest split primes > X
# genuine primes p_j ~ M^Cexp, = 1 mod 4, distinct, larger than X
lo = int(M ** Cexp); P = set()
while len(P) < r:
    c = random.randrange(lo, 2 * lo)
    if c % 4 == 1 and c > X and isprime(c): P.add(c)
P = sorted(P); pis = [gauss_prime(p) for p in P]
qs = list(primerange(3, X + 1))
designs = {k: [] for k in ('RAND', 'COH', 'TWIST', 'GEN')}
mods = []
for q in qs:
    T = Tq(q); m = T.m; mods.append(m)
    nu = np.array([1 if legendre_symbol(p % q, q) == -1 else 0 for p in P])
    # RAND
    lab = np.random.randint(0, m // 2, size=r) * 2 + nu
    designs['RAND'].append(lab % m)
    # COH
    coh = np.array([T.dlog_ratio(a, b) for (a, b) in zeta])
    designs['COH'].append(coh)
    # TWIST
    delta = (coh % 2) ^ nu
    if q % 8 in (3, 5):
        t = T.dlog_i()                      # i is a non-square iff q = +-3 mod 8
    else:
        ell = next(l for l in small if legendre_symbol(l % q, q) == -1)
        a, b = gauss_prime(ell); t = T.dlog_ratio(a, b)
    designs['TWIST'].append((coh + delta * t) % m)
    # GEN
    designs['GEN'].append(np.array([T.dlog_ratio(a, b) for (a, b) in pis]))
    assert all((designs['TWIST'][-1] % 2) == nu) and all((designs['GEN'][-1] % 2) == nu)
mods = np.array(mods); logq = np.log(np.array(qs, dtype=float)); quarter = mods // 4
logq = logq * (np.array(qs) >= M ** QMIN)
lM = math.log(M)
print(f"M={M} r={r} X={X} #moduli={len(qs)} p_j ~ M^{Cexp}  counting q >= M^{QMIN}")
for name, labs in designs.items():
    L = np.stack(labs)                      # (nq, r)
    res = {}
    def pu(V):
        """V: (nq, ...) values mod m. returns max over units s of per-unit demand, shape (...)."""
        Vm = V % mods.reshape((-1,) + (1,) * (V.ndim - 1))
        qq = quarter.reshape((-1,) + (1,) * (V.ndim - 1))
        hit = (Vm % qq == 0); cls = Vm // qq
        lq = logq.reshape((-1,) + (1,) * (V.ndim - 1))
        return np.max(np.stack([((hit & (cls == u)) * lq).sum(0) for u in range(4)]), axis=0)
    res[1] = pu(L).max() / lM
    b2 = pu(2 * L).max()
    for s_ in (1, -1):
        dem = pu(L[:, :, None] + s_ * L[:, None, :]); np.fill_diagonal(dem, 0); b2 = max(b2, dem.max())
    res[2] = b2 / lM
    b3 = 0.0
    if KMAX >= 3:
        for s1 in (1, -1):
            for s2 in (1, -1):
                for j in range(r):
                    dem = pu(L[:, j][:, None, None] + s1 * L[:, :, None] + s2 * L[:, None, :])
                    dem[j, :] = 0; dem[:, j] = 0; np.fill_diagonal(dem, 0); b3 = max(b3, dem.max())
    res[3] = b3 / lM
    print(f"  {name:5s}  max per-unit demand/log M:  k=1 {res[1]:6.2f}   k=2 {res[2]:6.2f}   k=3 {res[3]:6.2f}")
