"""(a) Lemma 5.2 closed form by exact enumeration (Fractions) for m <= 7, several rho.
(b) Arc quantisation: for a zero combination sum n_i sigma_i v_i = 0, sum n_i (pi - 2 d_i) is in (pi/2) Z.
    Check on all level-0 same-v certificates and ties of random sets (float, tolerance 1e-9)."""
import itertools, math, random
from fractions import Fraction
from rcore import rand_primes_like_author, gauss, Lprime, tight_list
bad = 0
for m in range(1, 8):
    for rho in (Fraction(1, 3), Fraction(1, 7), Fraction(2, 5)):
        tot = Fraction(0)
        for bits in itertools.product((0, 1), repeat=2 * m):
            lo, hi = bits[:m], bits[m:]
            if not any(lo[a] and hi[b] for a in range(m) for b in range(a + 1, m)):
                pr = Fraction(1)
                for x in bits: pr *= rho if x else 1 - rho
                tot += pr
        closed = (1 - rho) ** (m - 1) * (1 + (m - 1) * rho)
        bad += tot != closed
print('Lemma 5.2 closed form mismatches:', bad)
# quantisation
for (M, K, seed) in [(164, 3.5, 0), (164, 4.5, 1), (84, 3, 2)]:
    P = int(M ** K); ps = rand_primes_like_author(M - 1, P, 2 * P, random.Random(seed))
    gp = [gauss(p) for p in ps]; L = Lprime(2 * M, set(ps)); rels = tight_list(gp, L)
    byv = {}
    for v, s, sg, s2 in rels: byv.setdefault(v, []).append((s, sg, math.asin(math.sqrt(float(s2)))))
    arcs = []
    for v, lst in byv.items():
        for (s1, g1, d1), (s2, g2, d2) in itertools.combinations(lst, 2):
            if g1 != g2: arcs.append((2 * math.pi - 2 * d1 - 2 * d2) / (math.pi / 2))
    off = max((abs(a - round(a)) for a in arcs), default=0)
    print(f'M={M} K={K}: same-v opposite-orientation combos {len(arcs)}, arc/(pi/2) values {sorted(set(round(a) for a in arcs))}, max deviation {off:.1e}')
