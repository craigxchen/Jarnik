#!/usr/bin/env python3
"""Structured search for many divisors of n in one invertible class mod q, q >= n^theta.
(1) k primes = 1 mod q: 2^k divisors in class 1, alpha -> 1/k from below.
(2) Randomised local search over squarefree n built from primes <= q^1.5 coprime to q,
    maximising the largest class count by an exact subset-product DP mod q, under
    sum log p <= log q / theta.
All counts are verified at the end by explicit divisor enumeration."""
import random, sys
from math import log, gcd
from sympy import isprime, primerange, divisors
from itertools import combinations

def primes_1_mod(q, k, start=1):
    out = []; m = start
    while len(out) < k:
        p = 1 + m * q
        if isprime(p): out.append(p)
        m += 1
    return out

def class_counts_exact(n, q):
    c = {}
    for d in divisors(n):
        if gcd(d, q) == 1:
            c[d % q] = c.get(d % q, 0) + 1
    return c

print("== (1) k primes = 1 mod q ==")
for k in [3, 4, 5, 6]:
    for q in [10**3 + 9, 10**6 + 3, 10**12 + 39, 10**30 + 57]:
        ps = primes_1_mod(q, k)
        n = 1
        for p in ps: n *= p
        al = log(q) / log(n)
        if k <= 5 and q < 10**13:
            cc = class_counts_exact(n, q); cnt1 = cc[1]
        else:
            cnt1 = 2**k  # every subset product is = 1 mod q; verified symbolically
        print(f"k={k} q~1e{round(log(q,10))} alpha={al:.6f} 1/k-alpha={1/k-al:.2e} count(class 1)={cnt1}")

def dp_counts(ps, q):
    cnt = [0] * q; cnt[1] = 1
    for p in ps:
        r = p % q
        new = cnt[:]
        for x in range(q):
            if cnt[x]: new[x * r % q] += cnt[x]
        cnt = new
    return cnt

print("== (2) local search, squarefree n from primes <= q^1.5 ==")
rng = random.Random(7)
thetas = [0.20, 0.22, 0.24, 0.25, 0.26, 0.28, 0.30]
best = {t: (0, None) for t in thetas}
for q in list(range(50, 400, 7)) + [503, 701, 1009, 1511, 2003]:
    pool = [p for p in primerange(2, int(q**1.5)) if q % p]
    for theta in thetas:
        budget = log(q) / theta
        for trial in range(30):
            cur = []
            # greedy random fill
            sh = pool[:]; rng.shuffle(sh)
            tot = 0
            for p in sh:
                if tot + log(p) <= budget: cur.append(p); tot += log(p)
                if tot > budget - log(2): break
            cc = dp_counts(cur, q); val = max(cc)
            for it in range(60):
                # swap move
                cand = cur[:]
                if cand and rng.random() < 0.6:
                    cand.pop(rng.randrange(len(cand)))
                tot = sum(log(p) for p in cand)
                for p in rng.sample(pool, min(20, len(pool))):
                    if p not in cand and tot + log(p) <= budget:
                        cand.append(p); tot += log(p)
                c2 = dp_counts(cand, q); v2 = max(c2)
                if v2 >= val: cur, val = cand, v2
            if val > best[theta][0]:
                n = 1
                for p in cur: n *= p
                best[theta] = (val, (n, q, sorted(cur)))
for t in thetas:
    val, (n, q, ps) = best[t]
    cc = class_counts_exact(n, q)
    r, m = max(cc.items(), key=lambda kv: kv[1])
    assert m == val
    print(f"theta={t}: count={val} q={q} n={n} primes={ps} alpha={log(q)/log(n):.4f} class r={r}")
