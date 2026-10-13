# C7: the norm bound in Theorems A-C is necessary. For a column prime p and a bound X, a Gaussian prime
# zeta over a prime p' = p mod 4*prod{q <= X, q = 7 mod 8} is sub-genuine (p does not divide Norm zeta = p')
# yet its coherent label log(zeta/conj zeta) has the forced parity nu(q) = [(p/q) = -1] at EVERY prime
# q = 7 mod 8 below X: parity of the label is [(Norm zeta/q) = -1] (sharp.md Prop 2.1(3); re-tested
# exactly in t_glib.py) and (p'/q) = (p/q). So exact parity-coherence on Q_7 cap [1,X] with a sub-genuine
# zeta is possible once Norm zeta is allowed to be exp(O(X)).
from glib import *
import math
for p, X in [(13, 200), (101, 500), (1009, 1000)]:
    assert is_prime(p) and p % 4 == 1
    qs = [q for q in sieve(X) if q % 8 == 7]
    Pi = 4 * math.prod(qs)
    k = 1
    while not is_prime(p + k * Pi): k += 1
    pp = p + k * Pi
    ok = all(legendre(pp, q) == legendre(p, q) for q in qs)
    print(f"p={p}, X={X}: p' = p + {k}*Pi ({len(str(pp))} digits, prime, = 1 mod 4: {pp % 4 == 1}); "
          f"(p'/q) = (p/q) at all {len(qs)} primes q = 7 mod 8 <= X: {ok}")
