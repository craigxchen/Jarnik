import random
from fractions import Fraction as Fr
from check3 import marks, stats
random.seed(1)
for T in (1, 2):
    worst = 0; best = 1e9; cnt = 0
    for trial in range(60):
        seeds = set()
        while len(seeds) < 4:
            q = random.randint(1, 30); p = random.randint(q, 10*q)
            seeds.add(Fr(p, q))
        sd = [(s.numerator, s.denominator) for s in seeds]
        for n in (1, 3, 8, 20, 50):
            gn, H, DH, gam = stats(marks(T, n, sd)); cnt += 1
            worst = max(worst, DH); best = min(best, DH)
    print(f'T={T}: {cnt} tuples, seeds in [1,10] with denominators <=30: Delta*H^(1/4) in [{best:.4g}, {worst:.4g}]')
