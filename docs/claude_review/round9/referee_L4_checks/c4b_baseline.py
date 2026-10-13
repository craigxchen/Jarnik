# C4b: exact baseline for the per-column agreement notion of L4 Def. 3.1 (ell in u L + (m/4)Z, u a unit),
# against GENUINE labels G_j = log(pi_j/conj pi_j) on Q_7 cap [1, X].
# For x uniform in the parity class nu_j(q) of Z/m_Q and n = m_Q/4 (even), d = gcd(G_j(q), n):
#   P(agree) = #{y in Z/n of parity nu : gcd(y,n) = d} / (n/2) = phi(n/d) / (n/2)   (exact).
# We report the log-mass-weighted expectation per column (exact rational arithmetic per q, float at the end),
# i.e. the agreement that a uniformly random parity-respecting column has with the genuine column in expectation.
import math, random
from glib import *

def phi(n):
    r = n
    for l in factor(n): r = r // l * (l - 1)
    return r

X = 200000
q7 = [q for q in sieve(X) if q % 8 == 7]
tori = {}
for q in q7:
    A = 1
    while q ** (A + 1) <= X: A += 1
    tori[q] = Torus(q, A)
rng = random.Random(5)
cols = []
while len(cols) < 40:
    p = rng.randrange(10 ** 3, 10 ** 6)
    if p % 4 == 1 and is_prime(p) and p not in cols: cols.append(p)
tot = sum(math.log(q) for q in q7)
vals = []
for p in cols:
    pi = two_squares(p)
    e = 0.0
    for q in q7:
        T = tori[q]
        G = T.log(ratio(pi, T.Q))
        n = T.m // 4
        d = math.gcd(G, n)
        assert (G % 2) == (1 if legendre(p, q) == -1 else 0)
        e += math.log(q) * phi(n // d) / (n / 2)
    vals.append(e / tot)
vals.sort()
print(f"X={X}, {len(q7)} primes q=7 mod 8, 40 columns p in [1e3,1e6]: expected agreement of a uniformly random "
      f"parity-respecting column with the genuine column (Def. 3.1, unit u): min {vals[0]:.4f} median {vals[20]:.4f} max {vals[-1]:.4f}")
