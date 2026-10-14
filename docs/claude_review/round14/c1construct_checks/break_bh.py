# Thm B3: the systems of round10 L4K Theorem 2 (genuine below q_T, the vector (2.1) with a random permutation
# sigma_q above) are not fakes: at EVERY q > q_T, for EVERY sigma and parity pattern, many norm-4 monomials
# v = e_i + e_j - e_k - e_l have <v, a> = 0 exactly, hence collide with unit 0 at every level q^a <= X.
# K = 2, 3 (h = 1, beta_j = j, B = r).  Exact computation on an instance; explicit witness v.
import sys, math, random, itertools
sys.path.insert(0, '.')
from lib14 import split_primes, nu
from sympy import primerange
from collections import defaultdict

def avec(K, r, sig, nuq):
    B = r
    Lam = K * (K + 1) * B + K
    return [2 * (Lam + (K + 1) * sig[j]) + nuq[j] for j in range(r)]

def run(K, r, plo, phi, X, seed):
    rng = random.Random(seed)
    ps = sorted(rng.sample(split_primes(plo, phi), r))
    B = r
    AK = K * (2 * (K + 1)**2 * B + 2 * K + 1)
    qT = 4 * AK + 5
    quads = []  # unordered {i,j} vs {k,l}, disjoint, as (i,j,k,l) with i<j, k<l, (i,j)<(k,l)
    pairs = list(itertools.combinations(range(r), 2))
    for x in range(len(pairs)):
        for y in range(x + 1, len(pairs)):
            (i, j), (k, l) = pairs[x], pairs[y]
            if len({i, j, k, l}) == 4:
                quads.append((i, j, k, l))
    logm = defaultdict(float)   # log m_(v,0) restricted to q > qT
    hits = defaultdict(list)
    mincount = None
    nq = 0
    for q in primerange(qT + 1, X + 1):
        if q in ps: continue
        nq += 1
        sig = list(range(1, r + 1)); rng.shuffle(sig)   # beta_{sigma(j)} = sigma(j) in [1, r]
        nuq = [nu(p, q) for p in ps]
        a = avec(K, r, sig, nuq)
        Aq = int(math.floor(math.log(X) / math.log(q) + 1e-12))
        cnt = 0
        for (i, j, k, l) in quads:
            if a[i] + a[j] - a[k] - a[l] == 0:
                cnt += 1
                logm[(i, j, k, l)] += Aq * math.log(q)
                hits[(i, j, k, l)].append(q)
        mincount = cnt if mincount is None else min(mincount, cnt)
    # deterministic bound: >= (1/3) sum_s N_s^2 - ..., compare with r^3/9
    best = max(logm, key=lambda t: logm[t]) if logm else None
    lp = [math.log(p) for p in ps]
    budget = lambda t: 0.5 * (lp[t[0]] + lp[t[1]] + lp[t[2]] + lp[t[3]])
    viol = sum(1 for t in quads if logm[t] > budget(t))
    print(f"K={K} r={r} primes in [{plo},{phi}] X={X}: q_T={qT}, moduli q in (q_T,X]: {nq}, norm-4 monomials of the form e_i+e_j-e_k-e_l: {len(quads)}")
    Ns = {}
    for (x, y) in itertools.combinations(range(1, r + 1), 2):
        Ns[x + y] = Ns.get(x + y, 0) + 1
    Zr = sum(N * N for N in Ns.values()) / 6 - (r * (r - 1) // 2) / 2
    print(f"   proved lower bound Z_r = (1/6) sum_s N_s^2 - (1/2) binom(r,2) = {Zr:.2f}")
    print(f"   min over q of #such monomials with <v,a>=0 exactly: {mincount}  (total such monomials {len(quads)})")
    print(f"   monomials violating the unit-0 budget w(v)/2 from moduli > q_T alone: {viol} of {len(quads)}")
    if best:
        print(f"   witness v = e_{best[0]+1}+e_{best[1]+1}-e_{best[2]+1}-e_{best[3]+1} (primes {ps[best[0]]},{ps[best[1]]},{ps[best[2]]},{ps[best[3]]}):"
              f" log m_(v,0) >= {logm[best]:.1f} vs budget {budget(best):.2f}; collides at {len(hits[best])} moduli, e.g. {hits[best][:6]}")
    return ps

for (K, r, plo, phi, X, seed) in [(2, 10, 3000, 6000, 30000, 11), (3, 10, 3000, 6000, 60000, 12), (3, 16, 20000, 40000, 60000, 13)]:
    run(K, r, plo, phi, X, seed)
