#!/usr/bin/env python3
"""Independent check of stab.md Lemma 3.1 (cofactor) on non-simple profiles, and of the chain-deficit example.
F(zeta) = sum_m c_m zeta^{m - m0},  D_S = prod_p pi^{min<m,a>-<m0,a>} pibar^{-max<m,a>-<m0,e1-a>}.
Claims: D_S | F in Z[i];  N(F)/N(D_S) = |f(y)|^2 prod_p p^{width_S(a(p))}  (y_j = zeta_j/zeta_5);
width_S(a(p)) <= sum_t width_S(1_{S_t}) (subadditivity along the layer chain) -> chain deficit >= 0."""
import random, itertools, math, sys
import mpmath as mp
sys.path.insert(0, '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round27/kill_stab_checks')
from cox import gm, gc, gn, gpow, gdivexact, primes1mod4, gauss_prime, tuple_from_profile, random_profile, UNITS
mp.mp.dps = 60

def gpowz(x, k):
    return gpow(x, k)

def run(trials=300, seed=31):
    rnd = random.Random(seed)
    plist = primes1mod4(5, 200)
    fails_div = fails_norm = fails_sub = 0; maxdef = 0.0; n = 0
    for _ in range(trials):
        prof = random_profile(rnd, plist, rnd.randint(2, 5), rnd.choice([2, 3, 4]))
        eps = rnd.choice(UNITS)
        Z = tuple_from_profile(prof, eps)
        # random support in M (sum 0), small entries
        S = set()
        while len(S) < rnd.randint(2, 5):
            v = [rnd.randint(-2, 2) for _ in range(4)]; S.add(tuple(v + [-sum(v)]))
        S = sorted(S)
        coef = {m: (rnd.randint(-5, 5), rnd.randint(-5, 5)) for m in S}
        if all(c == (0, 0) for c in coef.values()): continue
        m0 = [min(m[j] for m in S) for j in range(5)]
        F = (0, 0)
        for m in S:
            term = coef[m]
            for j in range(5): term = gm(term, gpowz(Z[j], m[j] - m0[j]))
            F = (F[0] + term[0], F[1] + term[1])
        D = (1, 0); logW = 0.0; Wprod = 1
        for p, (a, e, pi) in prof.items():
            mn = min(sum(mi * ai for mi, ai in zip(m, a)) for m in S)
            mx = max(sum(mi * ai for mi, ai in zip(m, a)) for m in S)
            ea = -mx - sum(m0[j] * (e - a[j]) for j in range(5))
            pa = mn - sum(m0[j] * a[j] for j in range(5))
            assert pa >= 0 and ea >= 0
            D = gm(D, gm(gpow(pi, pa), gpow(gc(pi), ea)))
            Wprod *= p ** (mx - mn)
            # subadditivity along layer chain
            chain = 0
            for t in range(1, e + 1):
                St = [1 if a[j] >= t else 0 for j in range(5)]
                vals = [sum(mi * si for mi, si in zip(m, St)) for m in S]
                chain += max(vals) - min(vals)
            if mx - mn > chain: fails_sub += 1
            maxdef = max(maxdef, (chain - (mx - mn)) * math.log(p))
        n += 1
        c = gdivexact(F, D)
        if c is None: fails_div += 1; continue
        # |f(y)|^2 with y_j = zeta_j / zeta_5
        y = [mp.mpc(*Z[j]) / mp.mpc(*Z[4]) for j in range(5)]
        fy = mp.fsum(mp.mpc(*coef[m]) * mp.fprod(y[j] ** m[j] for j in range(5)) for m in S)
        lhs = mp.mpf(gn(c)); rhs = abs(fy) ** 2 * Wprod
        if lhs == 0 and abs(fy) < mp.mpf(10) ** -40: continue
        if abs(lhs - rhs) > mp.mpf(10) ** -30 * max(lhs, 1): fails_norm += 1
    print(f'Lemma 3.1: {n} cases (non-simple profiles, 2-4 layers/prime): divisibility fails {fails_div}, '
          f'norm identity fails {fails_norm}, subadditivity fails {fails_sub}; max per-prime chain deficit {maxdef:.2f}')

def example():
    m = (1, -1, 1, -1, 0); a = (2, 1, 0, 1, 0)
    layers = [[1 if a[j] >= t else 0 for j in range(5)] for t in (1, 2)]
    w_layers = [abs(sum(x * y for x, y in zip(m, L))) for L in layers]
    print('example a(p)=(2,1,0,1,0), m=e1-e2+e3-e4: layer sets', [[j + 1 for j in range(5) if L[j]] for L in layers],
          'layer widths', w_layers, 'width(a(p)) =', abs(sum(x * y for x, y in zip(m, a))),
          '-> deficit per prime = 2 log p')
    # max achievable deficit for this S over balanced profiles: 8 +-1 cuts pairable into nested opposite-sign chains
    CUTS = [frozenset(A) for r in range(1, 5) for A in itertools.combinations(range(1, 5), r)]
    val = {A: sum(m[j - 1] for j in A) for A in CUTS}
    print('cut values <m,1_A>:', {tuple(sorted(A)): v for A, v in val.items()}, ' sum|.| =', sum(abs(v) for v in val.values()))

if __name__ == '__main__':
    run(); example()
