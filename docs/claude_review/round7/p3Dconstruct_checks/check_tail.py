"""Evidence for Thm 1 / H1: single-column collision moduli log m^gen(e_j) at all odd prime powers <= X=3M,
for three designs on r = 8M columns with genuine parities nu_j(q) of random primes p_j ~ M^4:
  RAND : random fibre (labels uniform in parity class, random lifts),
  HYB  : coherent zeta_j + root repair on classes k(q)<=3, random on k(q)>=4,
  COH  : coherent zeta_j + root repair on all classes (exact Gaussian arithmetic).
Prints max_j and mean_j of log m^gen(e_j)/log M, and for COH the class split of the maximiser."""
import random, sys, math
from gauss import *

def run(M, rfac=8, seed=3):
    rng = random.Random(seed)
    X = 3*M; r = rfac*M
    ps = []
    while len(ps) < r:
        p = rng.randrange(M**4, 2*M**4)
        if p % 4 == 1 and is_prime(p): ps.append(p)
    hi = 4*X
    while True:
        zs = gaussian_primes_norm_range(X, hi)
        if len(zs) >= r: zs = zs[:r]; break
        hi *= 2
    qs = [q for q in primes_upto(X)[1:]]
    res = {d: [0.0]*r for d in ("RAND", "HYB", "COH")}
    split = [[0.0, 0.0] for _ in range(r)]       # COH: (k<=3 part, k>=4 part)
    for q in qs:
        A = 1
        while q**(A+1) <= X: A += 1
        k = v2(tq_order(q, 1)); lq = math.log(q)
        legp = [legendre(p, q) for p in ps]
        # random-style designs: level-1 hit prob from parity, then lifts with prob 1/q each
        for d in ("RAND", "HYB"):
            if d == "HYB" and k <= 3: continue
            m = tq_order(q, 1)
            for j in range(r):
                nu = 0 if legp[j] == 1 else 1
                lab = 2*rng.randrange(m//2) + nu
                if lab % (m//4) == 0:
                    a = 1
                    while a < A and rng.randrange(q) == 0: a += 1
                    res[d][j] += a*lq
        # coherent with root repair (exact), levels 1..A
        f, _ = sylow2_generator(q, A, rng)
        Q = q**A
        for j, (x, y) in enumerate(zs):
            eps = 0 if legendre(x*x+y*y, q) == legp[j] else 1
            z = (x % Q, y % Q)
            t = mul(mul(z, inv(conj(z, Q), Q), Q), pw(f, eps, Q), Q)
            a = 0
            for lev in range(1, A+1):
                Ql = q**lev
                tl = (t[0] % Ql, t[1] % Ql)
                if pw(tl, 4, Ql) == (1, 0): a = lev
                else: break
            if a:
                res["COH"][j] += a*lq
                split[j][0 if k <= 3 else 1] += a*lq
                if k <= 3: res["HYB"][j] += a*lq
    lM = math.log(M)
    out = {d: (max(v)/lM, sum(v)/len(v)/lM) for d, v in res.items()}
    jmax = max(range(r), key=lambda j: res["COH"][j])
    return out, (split[jmax][0]/lM, split[jmax][1]/lM)

if __name__ == "__main__":
    Ms = [int(a) for a in sys.argv[1:]] or [64, 128, 256, 512, 1024]
    print("M     r    | RAND max/mean      | HYB max/mean       | COH max/mean       | COH argmax split k<=3 / k>=4 | logM/loglogM")
    for M in Ms:
        out, sp = run(M)
        ll = math.log(M)/math.log(math.log(M))
        print(f"{M:5d} {8*M:6d} | {out['RAND'][0]:6.2f} {out['RAND'][1]:6.2f}      | {out['HYB'][0]:6.2f} {out['HYB'][1]:6.2f}      "
              f"| {out['COH'][0]:6.2f} {out['COH'][1]:6.2f}      | {sp[0]:5.2f} / {sp[1]:5.2f}                | {ll:5.2f}", flush=True)
