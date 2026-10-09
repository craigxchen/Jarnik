# Exhaustive check of the chain lemma:
# Phi(n) = sum_{tau=1}^{e} (S_tau - M/2)^2 + lam * sum_a n_a^2 >= c(lam) M^2/2,
# c(lam) = (sqrt(1+4 lam)-1)/2, lam = 1/(p-1); n ranges over compositions of M into e+1 parts.
from fractions import Fraction
from math import sqrt
from itertools import combinations

def compositions(M, parts):
    # stars and bars
    for bars in combinations(range(M + parts - 1), parts - 1):
        prev = -1; out = []
        for b in bars:
            out.append(b - prev - 1); prev = b
        out.append(M + parts - 1 - prev - 1)
        yield out

worst = 1e9
for p in [3, 5, 13, 29]:
    lam = Fraction(1, p - 1)
    c = (sqrt(1 + 4 * float(lam)) - 1) / 2
    for e in range(1, 6):
        for M in range(1, 19 if e <= 3 else 13):
            best = None
            for n in compositions(M, e + 1):
                S = 0; phi = Fraction(0)
                for tau in range(1, e + 1):
                    S += n[tau - 1]
                    phi += (Fraction(S) - Fraction(M, 2)) ** 2
                phi += lam * sum(x * x for x in n)
                if best is None or phi < best: best = phi
            ratio = float(best) / (c * M * M / 2)
            worst = min(worst, ratio)
            if ratio < 1 - 1e-12:
                print("VIOLATION", p, e, M, best)
    print("p", p, "done; worst ratio so far", worst)
