"""Referee (analytic lens): heavy-tail computation behind sharp.md Prop. 6.5, done per unit.

Random-fibre model (Lemma 6.3(2)): a non-twin column label ell_j^(q) is uniform on the m_q/2
residues of parity nu_j(q).  The single-prime monomial e_j collides at q with unit i^s iff
ell_j = iota s (mod m_q); iota s runs over the 4-torsion {0, m/4, m/2, 3m/4}.
  * aggregated (any unit):  Pr = 4/m_q                     (what Prop 6.5 counts)
  * one fixed unit s:       Pr = 2/m_q if parity(iota s) = nu_j(q), else 0
    -> averaged over a fair parity bit: 1/m_q, independent over q.
Height condition (Prop 6.1) is per unit:  log m_(e_j,s) <= w_j/2 + log(nu_s^-1) <= tau/2 + 0.35.

For X = 3M (level 1 only; prime powers add o(1)) we compute by exact DP over primes the law of
Y = sum_q log q * Bernoulli(p_q) (bins of 0.02, rounded UP so tails are not underestimated) and
the smallest tau with  E#{violating (column, unit)} = 4 r Pr[Y_unit > tau/2] <= 1  (per unit)
versus  r Pr[Y_agg > tau/2] <= 1  (aggregated, as in Prop 6.5), r = 33 M.  Units:
log^2 M/loglog M.  Prop 6.5 (corrected) predicts both tend to 2 (slowly).
Part 2: Monte Carlo of the actual max over r columns of the per-unit and aggregated demand.
"""
import math
import numpy as np

def primes_upto(n):
    s = np.ones(n + 1, dtype=bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return np.nonzero(s)[0]

h = 0.02
rng = np.random.default_rng(11)
print("Part 1 (exact DP tails):")
for M in (10 ** 3, 10 ** 4, 10 ** 5, 10 ** 6):
    r = 33 * M
    qs = [int(q) for q in primes_upto(3 * M) if q > 2]
    L = math.log(M)
    unit_scale = L * L / math.log(L)
    B = int(6 * unit_scale / h) + 10
    res = {}
    for label, factor, count in (("per-unit", 1.0, 4 * r), ("aggregated", 4.0, r)):
        dist = np.zeros(B); dist[0] = 1.0
        for q in qs:
            m = q - 1 if q % 4 == 1 else q + 1
            p = min(1.0, factor / m)
            k = int(math.ceil(math.log(q) / h))
            new = dist * (1 - p)
            new[k:] += p * dist[:-k]
            dist = new
        tail = np.cumsum(dist[::-1])[::-1]          # tail[i] = Pr[Y >= i h]
        idx = np.nonzero(count * tail <= 1.0)[0][0]  # smallest level with expected count <= 1
        tau = 2 * (idx * h - 0.35)
        res[label] = tau / unit_scale
    print(f"  M = {M:8d}: threshold tau / (log^2 M/loglog M):  per-unit {res['per-unit']:.3f},  aggregated {res['aggregated']:.3f}")

print("Part 2 (Monte Carlo max over r = 33 M columns, X = 3M, level 1):")
for M in (10 ** 3, 10 ** 4):
    r = 33 * M
    qs = [int(q) for q in primes_upto(3 * M) if q > 2]
    L = math.log(M)
    unit_scale = L * L / math.log(L)
    unit = np.zeros((4, r)); agg = np.zeros(r)
    for q in qs:
        m = q - 1 if q % 4 == 1 else q + 1
        nu = rng.integers(0, 2, size=r)                 # Legendre parity of p_j at q (fair bits)
        # 4-torsion elements s*m/4 and their parities
        par = [(s * (m // 4)) % 2 for s in range(4)]
        ell = 2 * rng.integers(0, m // 2, size=r) + nu  # uniform in parity class
        lq = math.log(q)
        anyhit = np.zeros(r, dtype=bool)
        for s in range(4):
            hit = ell == s * (m // 4)
            unit[s] += hit * lq
            anyhit |= hit
        agg += anyhit * lq
    nu_inv_log = np.array([0.0, 0.5 * math.log(2), 0.0, 0.5 * math.log(2)])[:, None]
    per = (unit - nu_inv_log).max()
    print(f"  M = {M:6d}: max_j log m_(e_j) [aggregated] = {agg.max() / unit_scale:.3f};"
          f"  max_(j,s) [log m_(e_j,s) - log nu_s^-1] = {per / unit_scale:.3f}  (units log^2 M/loglog M)")
