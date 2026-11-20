import pickle, random, math
from marks import gmul
ps, pis, q, marks = pickle.load(open("marks_200.pkl", "rb")); s = len(ps)
conj = lambda z: (z[0], -z[1])
f = lambda r, k: pis[k] if (r >> k) & 1 else conj(pis[k])
random.seed(5); bad = 0; n = 20000; minQsep = None
for _ in range(n):
    i, j = random.sample(range(len(marks)), 2)
    H = (1, 0); X = (1, 0); Q = 1
    for k in range(s):
        if ((i >> k) & 1) == ((j >> k) & 1): H = gmul(H, f(i, k))
        else: X = gmul(X, f(i, k)); Q *= ps[k]
    T = -X[1]
    d = (marks[j][0] - marks[i][0], marks[j][1] - marks[i][1])
    ok = d == gmul((0, 2), gmul(H, (T, 0))) and 0 < T * T < Q and math.gcd(*d) == 2 * abs(T) and X[0] != 0
    ok &= all(T % ps[k] for k in range(s) if ((i >> k) & 1) != ((j >> k) & 1))  # p does not divide T on separated pairs
    bad += not ok
print("random pairs checked:", n, "failures:", bad)
G = 4391525450
print("primes p | q with G^2 = 1 mod p (canonical squared-chord residue kept at p):", [p for p in ps if (G * G - 1) % p == 0])
print("C^4/16 =", 8**4/16, "<< q")
