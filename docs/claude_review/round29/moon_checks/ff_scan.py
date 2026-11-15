"""Function-field circle problem over F_q[t], q = 3 mod 4 (so i not in F_q, F_q[t][i] = F_{q^2}[t]).
Points: monic z in F_{q^2}[t] with z*zbar = N (bar = q-Frobenius on coefficients).
Ultrametric arc condition: deg(z_x - z_y) <= D/2 + c  <=>  equal coefficients in degrees > floor(D/2 + c).
For squarefree N = prod of distinct split primes P_j = pi_j pibar_j, monic points = 2^omega sign choices.
Max cluster at N = max fibre of z -> (coefficients of degree > D/2 + c).
Exhaustive over all such N with deg N <= nmax (for the given q).  Output: table of max M by (deg N, c).
"""
import sys, itertools
from collections import defaultdict

q = int(sys.argv[1]); nmax = int(sys.argv[2]); kmax = int(sys.argv[3])  # kmax = max degree of pi (P has degree 2k)
assert q % 4 == 3
Q2 = q * q
# element a + b*i  encoded as a + q*b
def enc(a, b): return (a % q) + q * (b % q)
def dec(x): return x % q, x // q
MUL = [[0] * Q2 for _ in range(Q2)]; ADD = [[0] * Q2 for _ in range(Q2)]
for x in range(Q2):
    a, b = dec(x)
    for y in range(Q2):
        c, d = dec(y)
        MUL[x][y] = enc(a * c - b * d, a * d + b * c)
        ADD[x][y] = enc(a + c, b + d)
NEG = [enc(-dec(x)[0], -dec(x)[1]) for x in range(Q2)]
CONJ = [enc(dec(x)[0], -dec(x)[1]) for x in range(Q2)]  # x^q = conj since i^q = -i

def pmul(f, g):
    r = [0] * (len(f) + len(g) - 1)
    for i, a in enumerate(f):
        if a == 0: continue
        for j, b in enumerate(g):
            if b: r[i + j] = ADD[r[i + j]][MUL[a][b]]
    return r
def pconj(f): return [CONJ[a] for a in f]
def is_fq(f): return all(dec(a)[1] == 0 for a in f)

# monic polys of degree k over F_{q^2}: coefficient list low->high, last = 1
def monics(k):
    for cs in itertools.product(range(Q2), repeat=k):
        yield list(cs) + [1]

# irreducible monic over F_{q^2} up to degree kmax (sieve by products)
irr = {}
for k in range(1, kmax + 1):
    red = set()
    # all products of lower-degree irreducibles of total degree k
    lower = [(d, f) for d in irr for f in irr[d]]
    def build(target, start, cur):
        if target == 0:
            red.add(tuple(cur)); return
        for idx in range(start, len(lower)):
            d, f = lower[idx]
            if d <= target:
                build(target - d, idx, pmul(cur, f))
    build(k, 0, [1])
    irr[k] = [f for f in monics(k) if tuple(f) not in red]
# split primes: pairs {pi, pibar} with pi != pibar (pi not over F_q)
split = []  # (degree of pi, pi)
seen = set()
for k in irr:
    for f in irr[k]:
        t = tuple(f)
        if t in seen: continue
        g = tuple(pconj(f))
        seen.add(t); seen.add(g)
        if g != t: split.append((k, f))
print(f"q={q}: split primes by deg(pi):", {k: sum(1 for d, _ in split if d == k) for k in irr}, flush=True)

best = defaultdict(int)  # (degN, c) -> max M
bestN = {}
cs = [-1, 0, 1]
count = 0
def process(chosen):
    global count
    count += 1
    pis = [split[i][1] for i in chosen]
    D = sum(split[i][0] for i in chosen)
    # enumerate points by gray-ish product
    pts = [[1]]
    for pi in pis:
        pb = pconj(pi)
        pts = [pmul(z, pi) for z in pts] + [pmul(z, pb) for z in pts]
    # remove the global conjugation duplicate? keep all (z and zbar are both points)
    for c in cs:
        thr = D / 2 + c  # keep coefficients of degree > thr
        lo = int(thr) + 1 if thr >= 0 else 0
        fib = defaultdict(int)
        for z in pts:
            fib[tuple(z[lo:D])] += 1
        m = max(fib.values())
        key = (2 * D, c)
        if m > best[key]:
            best[key] = m; bestN[key] = (tuple(chosen), D)

def rec(start, degleft, chosen):
    if chosen: process(chosen)
    for i in range(start, len(split)):
        d = split[i][0]
        if 2 * d <= degleft:
            chosen.append(i); rec(i + 1, degleft - 2 * d, chosen); chosen.pop()
rec(0, nmax, [])
print("N processed:", count)
for c in cs:
    row = [(n, best[(n, c)]) for n in range(2, nmax + 1, 2) if (n, c) in best]
    print(f"c={c:+d}: max M by deg N:", row)
# running max (since smaller degree clusters lift via composition only to multiples, report raw)
