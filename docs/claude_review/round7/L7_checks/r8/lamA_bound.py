# Rigorous upper bound for lambda_A(beta), exact over Z.
#   beta = (35/12) * avg_{sigma in S_7} sigma_*[gamma_A] + (1/12) * [mu]     (all coefficients >= 0)
#   lambda(C) <= comp(C) for every test curve C not contained in Sing(F_X)  (lamA_exact.lam_from_br)
# hence lambda_A(beta) <= (35/12) * avg_sigma comp(sigma(gamma_x)) + (1/12) * comp(mu).
# gamma_x = member of V_A through a random rational point x; its class is verified exactly (= gamma_A).
import sys, random, itertools, math, time
from fractions import Fraction as Fr
import sympy as sp
from star import member_through, analyze, point_of, LAB, t
from arr import PLANES, CHOSEN
import lamA_exact as E

GA = {('O','a1','a2'): 2, ('O','b1','b2'): 2, ('O','c1','c2'): 2, ('a1','b1','c1'): 1, ('a1','b1','c2'): 1,
      ('a1','b2','c1'): 1, ('a1','b2','c2'): 1, ('a2','b1','c1'): 1, ('a2','b2','c1'): 1, ('a2','b1'): 1,
      ('a2','b2'): 1, ('a2','c2'): 2, ('b1','c2'): 1, ('b2','c2'): 1}
pts = [point_of([PLANES[i], PLANES[j], PLANES[k]]) for (i, j, k) in CHOSEN]

def norm(mob):
    out = []
    for r in mob:
        g = 0
        for c in r: g = math.gcd(g, c)
        out.append([c // g for c in r])
    return out

def member_mob(seed):
    rng = random.Random(seed)
    x = [sp.Integer(c) for c in rng.sample(range(-6, 7), 7)]
    data, Q, _ = member_through(PLANES, pts, x)
    ev, inf, P = analyze(data)
    prod = sp.Poly(1, t, domain=sp.QQ)
    for e in ev.values(): prod = prod * e
    sqf = prod.sqf_part().degree() == prod.degree()
    cls = {tuple(LAB[a] for a in T): e.degree() for T, e in ev.items()}
    ok = Q.rank() == 4 and sqf and not inf and cls == GA
    mob = []
    for (N, D) in data:
        Np = sp.Poly(N, t, domain=sp.QQ); Dp = sp.Poly(D, t, domain=sp.QQ)
        co = [Np.coeff_monomial(t), Np.coeff_monomial(1), Dp.coeff_monomial(t), Dp.coeff_monomial(1)]
        den = 1
        for c in co: den = den * sp.Rational(c).q // math.gcd(den, sp.Rational(c).q)
        mob.append([int(c * den) for c in co])
    return ok, norm(mob), x

def coset_reps():
    idx = {l: a for a, l in enumerate(LAB)}
    cls = {frozenset(idx[l] for l in T): v for T, v in GA.items()}
    seen = {}
    for sig in itertools.permutations(range(7)):
        key = frozenset((frozenset(sig[a] for a in T), v) for T, v in cls.items())
        if key not in seen: seen[key] = sig
    return list(seen.values())

if __name__ == '__main__':
    seed0 = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    part, nparts = (int(sys.argv[2]), int(sys.argv[3])) if len(sys.argv) > 3 else (0, 1)
    random.seed(seed0)
    lm = E.lam_from_br(*E.leaf_brackets(E.int_mob(())))
    print('comp(mu) =', lm, flush=True)
    members = []
    s = seed0 * 1000
    while len(members) < 3:
        s += 1
        try:
            ok, mob, x = member_mob(s)
        except Exception as e:
            print('member construction failed for seed', s, e); continue
        print('member through x =', [str(c) for c in x], 'class verified:', ok, flush=True)
        if ok: members.append(mob)
    reps = coset_reps()
    print('number of classes sigma_*gamma_A:', len(reps), flush=True)
    mine = reps[part::nparts]
    tot = 0; hist = {}; fails = 0; t0 = time.time(); out = []
    for k, sig in enumerate(mine):
        val = None
        for mob in members:
            ms = [mob[sig.index(a)] for a in range(7)]     # label a of sigma(C) = label sigma^{-1}(a) of C
            try:
                val = E.lam_from_br(*E.leaf_brackets(ms)); break
            except Exception as e:
                fails += 1
        assert val is not None, ('no usable member for', sig)
        tot += val; hist[val] = hist.get(val, 0) + 1; out.append((sig, val))
        if k % 100 == 0: print(k, 'of', len(mine), 'hist', hist, 'fails', fails, 'time', round(time.time() - t0),
                               flush=True)
    print('PART', part, nparts, 'count', len(mine), 'sum', tot, 'hist', hist, 'fallbacks', fails, flush=True)
    with open('lamA_bound_part%d_of%d.txt' % (part, nparts), 'w') as f:
        f.write(repr({'comp_mu': lm, 'count': len(mine), 'sum': tot, 'hist': hist, 'vals': out}))
