"""Exact fraction of signed words of mass h (vectors v in Z^s, ||v||_1 = h) whose class-sum
vector is zero, for balanced partitions of s labels into b classes.  Referee's own code."""
import itertools
from fractions import Fraction


def vectors_of_mass(s, h):
    # compositions: choose support and signed parts
    for k in range(1, h + 1):
        for supp in itertools.combinations(range(s), k):
            for parts in itertools.product(range(1, h + 1), repeat=k):
                if sum(parts) != h:
                    continue
                for signs in itertools.product((1, -1), repeat=k):
                    yield {j: p * sg for j, p, sg in zip(supp, parts, signs)}


for s, b in [(24, 2), (24, 4), (24, 6), (48, 8)]:
    cls = [j % b for j in range(s)]
    out = []
    for h in (1, 2, 3):
        tot = bad = 0
        for v in vectors_of_mass(s, h):
            w = [0] * b
            for j, c in v.items():
                w[cls[j]] += c
            tot += 1
            bad += (not any(w))
        out.append((h, tot, bad, float(Fraction(bad, tot))))
    diff_pred = (s / b - 1) / (s - 1)
    print(f"s={s} b={b}: " + "; ".join(f"h={h} words={t} untested={u} frac={f:.4f}" for h, t, u, f in out)
          + f"; differences predicted (s/b-1)/(s-1)={diff_pred:.4f}")
