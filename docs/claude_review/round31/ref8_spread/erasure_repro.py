"""Independent reproduction of the 'COMPUTED (the erasure step)' numbers, and the K_C needed for C = 8.2394."""
import itertools, math
import numpy as np
from sympy import primerange

def K_of(q):
    sq = {x * x % q for x in range(1, q)}
    chi = np.array([0] + [1 if x in sq else -1 for x in range(1, q)])
    idx = (np.arange(q)[None, :] - np.arange(q)[:, None]) % q
    return chi[idx]   # K[i, j] = chi(j - i)

q = 127
K = K_of(q)
mins = {m: q for m in (2, 3, 4, 5, 7)}
words = []
for b in range(2, q):
    for d in range(2, q):
        if b != d:
            words.append({0: 1, b: 1, 1: -1, d: -1})   # a+b-c-d, a=0, c=1
for d in range(2, q):
    words.append({0: 2, 1: -1, d: -1})                # 2a-c-d, a=0, c=1
for w in words:
    v = np.zeros(q, dtype=np.int64)
    for k_, c in w.items():
        v[k_] += c
    Kv = K @ v
    for m in mins:
        mins[m] = min(mins[m], int(np.count_nonzero(Kv % m)))
print(f"q={q} primitive balanced length-4 words: min support of Kv mod m: {mins}")
v = np.zeros(q, dtype=np.int64); v[0] = 1; v[1] = -1
print(f"q={q} length-2 word e0-e1: support of Kv mod 2 = {np.count_nonzero((K @ v) % 2)}; isqrt(q-1)+1-h = {math.isqrt(q-1)+1-2}, sqrt(q)-h-1 = {q**.5-3:.2f}")

# Newton counts over F_11, t = 3: ordered pairs of disjoint 3-multisets, agreeing in first k power sums
p, t = 11, 3
ms = list(itertools.combinations_with_replacement(range(p), t))
ps = [tuple(sum(pow(x, s, p) for x in a) % p for s in range(1, t + 1)) for a in ms]
cnt = [0] * (t + 1)
for i, a in enumerate(ms):
    sa = set(a)
    for j, b in enumerate(ms):
        if sa & set(b):
            continue
        k = 0
        while k < t and ps[i][k] == ps[j][k]:
            k += 1
        for kk in range(k + 1):
            cnt[kk] += 1
print(f"F_{p}, t={t}: ordered disjoint pairs with first k moments equal, k=0..{t}: {cnt}")

# K_C for C = 8.2394: smallest initial segment of primes = 3 mod 4 with sum log l/(l+1) > log C
C = 8.23932686587
s = 0; S = []
for l in primerange(3, 10**7):
    if l % 4 == 3:
        s += math.log(l) / (l + 1); S.append(l)
        if s > math.log(C):
            break
print(f"C={C}: initial segment of primes 3 mod 4 up to {S[-1]} ({len(S)} primes) gives sum {s:.4f} > log C = {math.log(C):.4f}; K_C = {sum(S)}")
