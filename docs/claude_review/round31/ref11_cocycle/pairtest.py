"""Test of the congruence heuristic: for random odd v coprime to q in [V, 4V], count pairs of
eta=1 marks congruent mod v (equivalently v | T_ij), compare with C(N,2)/c(v),
c(v) = #{z mod v : |z|^2 = q mod v} = prod l^(e-1) (l - chi_{-4}(l))."""
import pickle, random, sys, math
from collections import Counter
from sympy import factorint
P = int(sys.argv[1]); V = int(float(sys.argv[2])); trials = int(sys.argv[3])
ps, pis, q, marks = pickle.load(open(f"marks_{P}.pkl", "rb"))
N = len(marks)
random.seed(1)
obs = 0; exp = 0.0
for _ in range(trials):
    while True:
        v = random.randrange(V, 4 * V) | 1
        if math.gcd(v, q) == 1:
            break
    c = 1
    for l, e in factorint(v).items():
        c *= l ** (e - 1) * (l - (1 if l % 4 == 1 else -1))
    cnt = Counter((x % v, y % v) for x, y in marks)
    obs += sum(n * (n - 1) // 2 for n in cnt.values())
    exp += N * (N - 1) / 2 / c
print(f"P={P} v in [{V},{4*V}) trials={trials}: congruent pairs observed {obs}, expected {exp:.1f}")
