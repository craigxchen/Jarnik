"""Integer lifting test (Data) on the actual (L_6)-refuting configurations Phi_N = Psi_t(Q1 + N(Q2-Q1)), t = t0.

A balanced 7-lift x = u/v (E_t frame) at scale w = w6/2 must satisfy u*beta_A - alpha_A*v = 0 mod e_A for every
contact A, with e_A | d_A (Sigma-free) and e_A ~ d_A^(1/2), and log max(|u|,|v|) = 8w + O(eta w + H).
For a split e = (e_A) the solutions form the lattice Lambda(e) = <(M,0),(r,1)>, M = prod e_A.
We record ell(e) = log lambda_1(Lambda(e)) (Euclidean, Gauss-reduced) in units of w, over sampled balanced splits.
A balanced bounded-residue 7-lift (hence any 8-lift, Lemma 1.1) needs ell(e)/w <= 8 + o(1) for some e.
Sanity check: each mover m of Phi_N itself satisfies the congruences with the FULL d_A (A ∋ m)."""
import sys, json, random, time
from fractions import Fraction as Fr
from math import log, gcd
import flint
from ell import quotient_curve
from lift import lift_T
from profile import profile

SIGMA = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 37, 41, 43, 61, 67, 71, 73, 97, 113, 311, 521, 571, 599, 1321, 2281]

def sigma_free(n):
    for p in SIGMA:
        while n % p == 0:
            n //= p
    return n

def gauss_reduce(b1, b2):
    def dot(x, y): return x[0] * y[0] + x[1] * y[1]
    if dot(b1, b1) > dot(b2, b2): b1, b2 = b2, b1
    while True:
        mu = Fr(dot(b1, b2), dot(b1, b1))
        m = round(mu)
        b2 = (b2[0] - m * b1[0], b2[1] - m * b1[1])
        if dot(b2, b2) >= dot(b1, b1):
            return b1, b2
        b1, b2 = b2, b1

def crt_r(conds):
    """conds: list of (modulus e, value x=alpha/beta as Fraction) -> (M, r) with r = x mod e for all."""
    M, r = 1, 0
    for e, x in conds:
        if e == 1: continue
        xr = (x.numerator % e) * pow(x.denominator % e, -1, e) % e
        # combine r mod M with xr mod e
        t = ((xr - r) * pow(M % e, -1, e)) % e
        r = r + M * t; M *= e
        r %= M
    return M, r

def divisors_from_fac(fac):
    ds = [1]
    for p, k in fac:
        ds = [d * p ** i for d in ds for i in range(k + 1)]
    return ds

def run(Nmax, nsample=4000, eta_s=0.25, out='intlift.json'):
    xs = [Fr(2), Fr(3), Fr(-5), Fr(1)]
    F, C, g, pts = quotient_curve(xs); k = g[3]
    Q1, Q2 = pts[0], pts[1]
    Dl = C.add(Q2, (Q1[0], -Q1[1])); R = Q1
    results = []
    for N in range(1, Nmax + 1):
        R = C.add(R, Dl)
        L = lift_T(F, R[0] / k)
        if not L:
            continue
        pt = L[0]
        try:
            d, b, Dm = profile(xs, pt)
        except ZeroDivisionError:
            continue
        sides = list(d)
        val = {Cs: (xs[min(Cs)] if min(Cs) <= 2 else xs[3]) for Cs in sides}
        dfree = {Cs: sigma_free(d[Cs]) for Cs in sides}
        logs = [log(x) for x in dfree.values() if x > 1]
        w6 = sum(log(x) for x in dfree.values()) / 25
        w = w6 / 2
        # sanity: movers satisfy full congruences
        san = True
        for mi, z in zip((3, 4, 5), pt):
            z = Fr(z)
            for Cs in sides:
                if mi in Cs and dfree[Cs] > 1:
                    e = dfree[Cs]; x = val[Cs]
                    if (z.numerator * x.denominator - x.numerator * z.denominator) % e != 0:
                        san = False
        # structure-beats-count for mover a (index 3)
        za = Fr(pt[0])
        Ma, ra = crt_r([(dfree[Cs], val[Cs]) for Cs in sides if 3 in Cs])
        b1, b2 = gauss_reduce((Ma, 0), (ra, 1))
        lam_a = log(max(abs(b1[0]), abs(b1[1])))
        h_a = log(max(abs(za.numerator), abs(za.denominator)))
        # factor
        t0 = time.time()
        facs = {Cs: [(int(p), int(e)) for p, e in flint.fmpz(dfree[Cs]).factor()] if dfree[Cs] > 1 else [] for Cs in sides}
        tf = time.time() - t0
        bal = {}
        for Cs in sides:
            D = dfree[Cs]; lD = log(D) if D > 1 else 0.0
            bal[Cs] = [e for e in divisors_from_fac(facs[Cs]) if D > 1 and abs(log(e) / lD - 0.5) <= eta_s]
        nbal = {Cs: len(v) for Cs, v in bal.items()}
        nonsplittable = sum(1 for v in nbal.values() if v == 0)
        rec = dict(N=N, w6=w6, w=w, sanity_movers=san, mover_a=dict(h_a_over_w=h_a / w, lambda1_over_w=lam_a / w,
                   index_over_w=log(Ma) / w), factor_time=tf, omega=[len(facs[Cs]) for Cs in sides],
                   nonsplittable_sides=nonsplittable)
        if True:
            # unsplittable sides are dropped (e_A = 1): this only weakens the necessary condition
            for Cs in sides:
                if not bal[Cs]:
                    bal[Cs] = [1]
            rnd = random.Random(N)
            ells = []
            skews = []
            total = 1
            for v in bal.values(): total *= len(v)
            if total <= nsample:
                import itertools
                combos = itertools.product(*[bal[Cs] for Cs in sides])
            else:
                combos = ([rnd.choice(bal[Cs]) for Cs in sides] for _ in range(nsample))
            for es in combos:
                M, r = crt_r([(e, val[Cs]) for e, Cs in zip(es, sides)])
                b1, b2 = gauss_reduce((M, 0), (r, 1))
                ells.append((log(max(abs(b1[0]), abs(b1[1]))) / w, log(M) / w))
                skews.append((log(M) / 2 - log(max(abs(b1[0]), abs(b1[1])))) / w)
            rec.update(n_sides_used=25 - nonsplittable, n_splits_total=total, n_tested=len(ells), min_ell_over_w=min(e for e, _ in ells),
                       mean_ell_over_w=sum(e for e, _ in ells) / len(ells),
                       mean_index_over_w=sum(i for _, i in ells) / len(ells),
                       max_skew_over_w=max(skews), mean_skew_over_w=sum(skews) / len(skews),
                       mover_a_skew_over_w=(log(Ma) / 2 - lam_a) / w)
        results.append(rec)
        print(json.dumps(rec), flush=True)
        json.dump(results, open(out, 'w'), indent=1)
    return results

if __name__ == '__main__':
    run(int(sys.argv[1]) if len(sys.argv) > 1 else 6)
