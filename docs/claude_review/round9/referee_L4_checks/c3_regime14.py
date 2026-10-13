# C3/C4: Construction D in a regime-(1.4)-type instance (X >> p_j, zeta_j = pi_sigma(j)), versus genuine labels.
# (a) per-unit height condition m_(v,s) <= nu_s^(-1) e^(w(v)/2) for all v of norm <= 2 (exact integers);
# (b) per-column agreement (Def. 3.1 sense: ell_j in u L + (m/4)Z, some unit u) with the GENUINE labels on Q_7.
import math, random, sys
from glib import *
from consD import *

def gcd_agree(x, y, n):
    return math.gcd(x % n, n) == math.gcd(y % n, n)

def run(X, r, seed, pmin, pmax):
    rng = random.Random(seed)
    cand = [p for p in range(pmin, pmax) if p % 4 == 1 and is_prime(p) and p != 17]
    ps = sorted(rng.sample(cand, r))
    pat = lambda p: tuple(legendre(p, q) for q in S0)
    classes = {}
    for j, p in enumerate(ps): classes.setdefault(pat(p), []).append(j)
    sigma = list(range(r))
    for cl in classes.values():
        for t, j in enumerate(cl): sigma[j] = cl[(t + 1) % len(cl)]
    moved = sum(1 for j in range(r) if sigma[j] != j)
    pis = [two_squares(p) for p in ps]
    S = System(ps, X)
    G = S.coherent(pis)
    systems = {
        'genuine': (G, pis),
        'D(id)': (S.construct_D(pis)[0], pis),
        'D(sigma)': (S.construct_D([pis[sigma[j]] for j in range(r)])[0], [pis[sigma[j]] for j in range(r)]),
    }
    print(f"X={X} r={r} seed={seed} primes={ps} ; sigma moves {moved}/{r} columns")
    mons = monomials(r, 1) + monomials(r, 2)
    for name, (ell, zs) in systems.items():
        viol = {1: 0, 2: 0}; worst = 0.0; nvs = {1: 0, 2: 0}; viol_v = {1: 0, 2: 0}
        for v in mons:
            k = sum(abs(c) for c in v.values())
            P = 1
            for j, c in v.items(): P *= ps[j] ** abs(c)
            prof = collision_profile(S, ell, v)
            anyv = False
            for s in range(4):
                m = mod_from(prof[s])
                nvs[k] += 1
                bound2 = P if s % 2 == 0 else 2 * P
                if m * m > bound2:
                    viol[k] += 1; anyv = True
                budget = 0.5 * math.log(P) + (0 if s % 2 == 0 else 0.5 * math.log(2))
                worst = max(worst, math.log(m) / budget)
            if anyv: viol_v[k] += 1
        print(f"  {name:9s}: height violations (v,s): norm1 {viol[1]}/{nvs[1]}, norm2 {viol[2]}/{nvs[2]};"
              f" monomials violating: norm1 {viol_v[1]}/{r}, norm2 {viol_v[2]}/{len(mons)-r};"
              f" max log m_(v,s)/budget = {worst:.2f}")
    # agreement with genuine labels on Q_7, log-mass weighted
    q7 = [q for q in S.qs if q % 8 == 7]
    tot = sum(math.log(q) for q in q7)
    for name in ('D(id)', 'D(sigma)'):
        ell = systems[name][0]
        fr = []
        for j in range(r):
            ag = sum(math.log(q) for q in q7 if gcd_agree(ell[q][j], G[q][j], S.T[q].m // 4))
            fr.append(ag / tot)
        print(f"  agreement of {name} with genuine G on Q_7 (log-mass): min {min(fr):.3f} median {sorted(fr)[r//2]:.3f} max {max(fr):.3f}")
    # control: random parity-respecting labels vs G
    fr = []
    for j in range(r):
        ag = 0.0
        for q in q7:
            T = S.T[q]
            x = rng.randrange(T.m // 2) * 2 + S.nu[q][j]
            if gcd_agree(x, G[q][j], T.m // 4): ag += math.log(q)
        fr.append(ag / tot)
    print(f"  agreement of random parity-respecting labels with G on Q_7: median {sorted(fr)[r//2]:.3f}")

for X, r, seed, pmin, pmax in [(20000, 12, 11, 29, 400), (50000, 10, 12, 29, 300)]:
    run(X, r, seed, pmin, pmax)
    sys.stdout.flush()
