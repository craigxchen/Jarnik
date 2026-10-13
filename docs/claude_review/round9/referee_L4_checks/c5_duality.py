# C5: exact check of the duality step of Lemma 4.2 as an F_2 linear system (any vector of
# quadratic characters at the odd primes l | K is realised by CRT), on random instances including
# ones with multiplicative relations; then evidence: an actual prime q = 7 mod 8 with all (n_c/q) = -1.
import random, math
from glib import *

def sqfree_odd_part(n):
    f = factor(n); m = 1
    for l, e in f.items():
        if l != 2 and e % 2 == 1: m *= l
    return m

def solvable(rows, rhs, nvars):
    # Gaussian elimination over F_2; rows are int bitmasks
    rows = list(rows); rhs = list(rhs)
    piv_row = 0
    for col in range(nvars):
        sel = next((i for i in range(piv_row, len(rows)) if rows[i] >> col & 1), None)
        if sel is None: continue
        rows[piv_row], rows[sel] = rows[sel], rows[piv_row]
        rhs[piv_row], rhs[sel] = rhs[sel], rhs[piv_row]
        for i in range(len(rows)):
            if i != piv_row and rows[i] >> col & 1:
                rows[i] ^= rows[piv_row]; rhs[i] ^= rhs[piv_row]
        piv_row += 1
    return all(not (rows[i] == 0 and rhs[i]) for i in range(len(rows)))

rng = random.Random(7)
gnorms = [a * a + b * b for a in range(0, 30) for b in range(0, 30) if a * a + b * b > 0]
qs7 = [q for q in sieve(2 * 10 ** 6) if q % 8 == 7]
fails, rel_inst, worstq, notfound = 0, 0, 0, 0
N_INST = 3000
for trial in range(N_INST):
    p = rng.choice([5, 13, 29, 37, 41, 53, 61, 73])
    P = rng.randrange(1, 9)
    norms = [rng.choice([g for g in gnorms if g % p != 0]) for _ in range(P)]
    ns = [p * a for a in norms]
    ms = [sqfree_odd_part(n) for n in ns]
    primes = sorted(set(l for m in ms for l in factor(m)))
    idx = {l: i for i, l in enumerate(primes)}
    rows = [sum(1 << idx[l] for l in factor(m)) for m in ms]
    rhs = [1 if -((-1) ** ((m - 1) // 2)) == -1 else 0 for m in ms]   # t_c = -sigma_(m_c); bit = [t_c = -1]
    # relation present?
    has_rel = any(math.isqrt(math.prod(ms[c] for c in range(P) if mask >> c & 1)) ** 2 ==
                  math.prod(ms[c] for c in range(P) if mask >> c & 1) for mask in range(1, 2 ** P))
    rel_inst += has_rel
    if not solvable(rows, rhs, len(primes)): fails += 1
    # evidence: least actual prime q = 7 mod 8 not dividing prod n_c with (n_c/q) = -1 for all c
    q = next((q for q in qs7 if all(n % q for n in ns) and all(jacobi(n % q, q) == -1 for n in ns)), None)
    if q is None: notfound += 1
    else: worstq = max(worstq, q)
print(f"Lemma 4.2 duality (F_2 system, exact): {N_INST} instances, {rel_inst} with nontrivial relations, "
      f"unsolvable: {fails}")
print(f"  evidence: witness prime q found below 2e6 in {N_INST - notfound} instances; largest least witness {worstq}")

worst = 0
for trial in range(3000):
    p = rng.randrange(10 ** 5, 10 ** 7)
    while not (p % 4 == 1 and is_prime(p)): p += 1
    a, b = rng.randrange(0, 60), rng.randrange(1, 60)
    n = p * (a * a + b * b)
    q = next(q for q in qs7 if n % q and jacobi(n % q, q) == -1)
    worst = max(worst, q)
print(f"Lemma 2.1 evidence: 3000 random sub-genuine n = p Norm(zeta), least q=7 mod 8 with (n/q)=-1 is at most {worst}")
