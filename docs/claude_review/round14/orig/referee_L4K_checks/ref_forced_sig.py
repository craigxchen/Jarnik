"""Independent exhaustive count of good level-1 systems on the forced set, K=2, ps={5,13,17,29}.
Signature of l in (Z/m)^r (parity-respecting) = set of (v,s), v in V_2 mod sign, with <v,l> = s*m/4 mod m.
A tuple of signatures is good iff for every (v,s): (prod of q with (v,s) in sig_q)^2 <= Norm(Z_v) * (2 if s odd).
Also: number of tuples with the genuine signature everywhere, and per-q genuine-signature vectors."""
import itertools, sys, time
import numpy as np
from collections import Counter
from rlib import *

ps = [5, 13, 17, 29]; K = 2; r = 4
F = [3, 7, 11, 19, 23, 31, 41]
V = VK_mod_sign(r, K)
pis = gaussian_primes(ps)
keys = [(i, s) for i in range(len(V)) for s in range(4)]
kidx = {k: t for t, k in enumerate(keys)}
bound2 = []
for (i, s) in keys:
    z = Zv(pis, V[i]); nz = z[0] ** 2 + z[1] ** 2
    bound2.append(nz * (2 if s % 2 else 1))
Varr = np.array(V)
per_q = []
gen_sig = {}
for q in F:
    m = q - chi4(q); n = m // 4
    nv = nu(ps, q)
    T = Torus(q)
    io = T.iota()
    L = [T.ratio_log(pi) for pi in pis]
    if io != m // 4:  # normalise so that iota = m/4
        L = [(-x) % m for x in L]
    choices = [[x for x in range(m) if x % 2 == nv[j]] for j in range(r)]
    allv = np.array(list(itertools.product(*choices)), dtype=np.int64)  # (count, r)
    t = (allv @ Varr.T) % m  # (count, 20)
    sigs = Counter()
    # encode signature as bitmask over 80 keys
    coll = (t % n == 0)
    sval = (t // n) % 4
    codes = np.zeros(len(allv), dtype=object)
    enc = [0] * len(allv)
    for i in range(len(V)):
        idx = np.nonzero(coll[:, i])[0]
        for row in idx:
            enc[row] |= 1 << kidx[(i, int(sval[row, i]))]
    sigs = Counter(enc)
    # genuine signature
    tg = [sum(a * b for a, b in zip(v, L)) % m for v in V]
    g = 0
    for i, tv in enumerate(tg):
        if tv % n == 0: g |= 1 << kidx[(i, (tv // n) % 4)]
    gen_sig[q] = g
    # cross-check genuine signature with D_s
    g2 = 0
    for i, v in enumerate(V):
        z = Zv(pis, v)
        for s in range(4):
            if Ds(z, s) % q == 0: g2 |= 1 << kidx[(i, s)]
    assert g == g2, q
    per_q.append((q, list(sigs.items())))
    print(f"q={q} parity vectors={len(allv)} signatures={len(sigs)} genuine-signature vectors={sigs[g]}", flush=True)

# DFS, larger q first
order = sorted(range(len(F)), key=lambda i: -F[i])
prods = [1] * len(keys)
count = 0; ntup = 0; gen_count = 0
nongen_at = Counter()
t0 = time.time()
def bits(x):
    out = []
    while x:
        b = x & -x; out.append(b.bit_length() - 1); x ^= b
    return out
sig_bits = {q: [(code, mult, bits(code)) for code, mult in lst] for q, lst in per_q}
def dfs(d, mult, allgen, nongen):
    global count, ntup, gen_count
    if d == len(order):
        count += mult; ntup += 1
        if allgen: gen_count += mult
        for q in nongen: nongen_at[q] += mult
        return
    q = F[order[d]]
    for code, mu, bl in sig_bits[q]:
        ok = True
        for k in bl:
            pq = prods[k] * q
            if pq * pq > bound2[k]: ok = False; break
        if not ok: continue
        for k in bl: prods[k] *= q
        isg = (code == gen_sig[q])
        dfs(d + 1, mult * mu, allgen and isg, nongen if isg else nongen + (q,))
        for k in bl: prods[k] //= q
dfs(0, 1, True, ())
print("good level-1 systems:", count, " signature tuples:", ntup)
print("all-genuine-signature systems:", gen_count, " = 2^%.3f" % (math.log2(gen_count) if gen_count else -1), " fraction %.3e" % (gen_count / count))
for q in F:
    print(f"  share of good systems with non-genuine signature at q={q}: {nongen_at[q]/count:.3f}")
print("time %.1fs" % (time.time() - t0))
