# Evidence (no proof uses it): demand of small monomials on the DEEP classes k(q) > E0, level 1, q <= X,
# for three designs: (R) coherent + root repair (zeta-norms = 1 mod 2^(E0+1)), (U) independent uniform
# parity-respecting labels, (G) genuine labels.  Reports mean and max over monomials of norm 1 and 2 of the
# deep demand  sum_{q deep, v collides at q (some unit)} log q,  and the certified bound sum_k log|Im gamma^(2^k)|
# restricted to k > E0 (which grows like 2^k).
import sys, math, random, itertools
sys.path.insert(0, '.')
from lib14 import *
from sympy import primerange
import lib10

def two_sylow_gen(q, k, rng):
    m = mq(q)
    while True:
        a, b = rng.randrange(q), rng.randrange(q)
        if (a * a + b * b) % q == 0: continue
        t = tor((a, b), q)
        f = mpow(t, m // (2 ** k), q)
        if mpow(f, 2 ** (k - 1), q) != (1, 0):
            return f

def run(r, E0, X, seed):
    rng = random.Random(seed)
    ps = sorted(rng.sample(split_primes(10**4, 3 * 10**4), r))
    ls = [l for l in primerange(3, 10**6) if l % (2 ** (E0 + 1)) == 1 and l not in ps][:r]
    zs = [gauss_prime(l) for l in ls]
    pis = [gauss_prime(p) for p in ps]
    vs = lib10.monos(r, 2, half=True)
    dem = {d: {(v, u): 0.0 for v in vs for u in range(4)} for d in "RUG"}
    nq = 0
    for q in primerange(3, X + 1):
        k = kq(q)
        if k <= E0 or q in ps or q in ls: continue
        nq += 1; m = mq(q); n4 = m // 4
        f = two_sylow_gen(q, k, rng)
        tR = []; tG = []
        for j in range(r):
            eps = 1 if legendre_symbol(ls[j] % q, q) != legendre_symbol(ps[j] % q, q) else 0
            t = tor(zs[j], q)
            if eps: t = mmul(t, f, q)
            tR.append(t); tG.append(tor(pis[j], q))
        # uniform parity-respecting labels: abstract Z/m labels, collision iff <v,ell> = 0 mod n
        nus = [nu(p, q) for p in ps]
        lab = [2 * rng.randrange(m // 2) + nus[j] for j in range(r)]
        units = [unit_mod(s, q) for s in range(4)]
        for v in vs:
            for d, ts in (("R", tR), ("G", tG)):
                tv = (1, 0)
                for j in range(r):
                    if v[j]:
                        tv = mmul(tv, mpow(ts[j] if v[j] > 0 else minv(ts[j], q), abs(v[j]), q), q)
                if tv in units: dem[d][(v, units.index(tv))] += math.log(q)
            tot = sum(v[j] * lab[j] for j in range(r)) % m
            if tot % n4 == 0: dem["U"][(v, tot // n4)] += math.log(q)
    print(f"r={r} E0={E0} X={X}: deep moduli (k(q) > E0, level 1): {nq}; zeta-norms {ls[:4]}...; p_j in [{ps[0]},{ps[-1]}]")
    for normk in (1, 2):
        sel = [v for v in vs if sum(map(abs, v)) == normk]
        wv_ = lambda v: sum(abs(t) * math.log(p) for t, p in zip(v, ps))
        line = f"   norm {normk} ({len(sel)} monomials x 4 units; per-unit deep demand / per-unit budget):"
        for d, name in (("R", "coherent+repair"), ("U", "uniform random"), ("G", "genuine")):
            rat = [dem[d][(v, u)] / (wv_(v) / 2 + (0.5 * math.log(2) if u % 2 else 0)) for v in sel for u in range(4)]
            over = sum(1 for x in rat if x > 1)
            line += f"  {name}: mean {sum(rat)/len(rat):.3f} max {max(rat):.2f} #>1: {over};"
        print(line)
    # certified bound for the deep classes of the coherent design (k from E0+1 to max k of q <= X)
    kmax = max(kq(q) for q in primerange(3, X + 1))
    print(f"   certified bound for coherent+repair on deep classes ~ sum_(k=E0+1)^(kmax={kmax}) 2^(k-1) log l per unit norm = {sum(2**(k-1) for k in range(E0+1, kmax+1))} * log l")

run(10, 2, 200000, 1)
run(10, 3, 200000, 2)
