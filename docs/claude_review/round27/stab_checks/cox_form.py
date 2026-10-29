"""Exact checks of the elementary Gaussian Cox form of 5-tuples (round 27, lane stab).

For random genuine Gaussian profiles (primes p = 1 mod 4, exponent profiles a_j(p) with min 0, max e_p):
 (1) build zeta_j = eps prod_p pi^{a_j} pibar^{e-a_j}; check primitivity (Gaussian gcd 1);
 (2) build gamma_A (A subset [4]) from threshold layers; check zeta_a = eps prod_{A ni a} g_A prod_{A not ni a} conj g_A,
     zeta_5 = eps prod conj g_A; chain condition (Ch);
 (3) Pi_a = prod_{A ni a} g_A: zeta_a conj(Pi_a) = zeta_5 Pi_a; |zeta_a - zeta_5|^2 = 4 R''^2 (Im Pi_a)^2/N(Pi_a);
 (4) pair cofactors Q_ab, t_ab = Im Q_ab, |zeta_a-zeta_b|^2 = 4 R''^2 t_ab^2 e^{-d_ab}; c_ab = (zeta_a-zeta_b)/gcd = unit*2i t_ab;
 (5) triangle relations Q_ab Q_bc = G^(b)_abc Q_ac for all ordered triples of [5];
 (6) Im/Re forms x_b t_a - x_a t_b = G_ab t_ab, x_a x_b + t_a t_b = G_ab y_ab, quadratic-form identity;
 (7) converse: random gamma's violating (Ch) give imprimitive tuples; random gamma's satisfying (Ch) give primitive ones;
 (8) unit classes: mixing units breaks integrality of t_ab (|c_ab|^2 = 2 mod 4 cases), as in kill_gauss minor 2.
"""
import random, sys
from itertools import combinations, permutations
from gauss_lib import *

FIVE = [1, 2, 3, 4, 5]
CUTS = [frozenset(s) for k in range(1, 5) for s in combinations([1, 2, 3, 4], k)]  # 15 sides not containing 5
assert len(CUTS) == 15

def build_tuple(primes, rng, emax=3):
    eps = (1, 0)
    prof = {}
    for p in primes:
        e = rng.randint(1, emax)
        a = [rng.randint(0, e) for _ in range(5)]
        i0, i1 = rng.sample(range(5), 2)
        a[i0] = 0; a[i1] = e
        prof[p] = (e, a)
    zeta = []
    for j in range(5):
        z = eps
        for p in primes:
            pi = gauss_prime(p)
            e, a = prof[p]
            z = gmul(z, gmul(gpow(pi, a[j]), gpow(gconj(pi), e - a[j])))
        zeta.append(z)
    return eps, prof, zeta

def gammas_from_profile(primes, prof):
    g = {A: (1, 0) for A in CUTS}
    chains = {}
    for p in primes:
        pi = gauss_prime(p)
        e, a = prof[p]
        chains[p] = ([], [])
        for t in range(1, e+1):
            S = frozenset(j+1 for j in range(5) if a[j] >= t)
            if 5 not in S:
                g[S] = gmul(g[S], pi); chains[p][0].append(S)
            else:
                A = frozenset(FIVE) - S
                g[A] = gmul(g[A], gconj(pi)); chains[p][1].append(A)
    return g, chains

def chain_ok(ch):
    pis, pibs = ch
    def is_chain(L):
        for X, Y in combinations(L, 2):
            if not (X <= Y or Y <= X):
                return False
        return True
    U1 = set().union(*pis) if pis else set()
    U2 = set().union(*pibs) if pibs else set()
    return is_chain(pis) and is_chain(pibs) and not (U1 & U2)

def zeta_from_gammas(eps, g):
    zs = []
    for a in [1, 2, 3, 4]:
        z = eps
        for A in CUTS:
            z = gmul(z, g[A] if a in A else gconj(g[A]))
        zs.append(z)
    z5 = eps
    for A in CUTS:
        z5 = gmul(z5, gconj(g[A]))
    zs.append(z5)
    return zs

def side(A, a):
    return a in A  # a in [5]; 5 never in A

def Qpair(g, a, b):
    """pair cofactor: prod_{A: a in A, b not} g_A * prod_{A: b in A, a not} conj g_A"""
    q = (1, 0)
    for A in CUTS:
        ia, ib = (a in A), (b in A)
        if ia and not ib:
            q = gmul(q, g[A])
        elif ib and not ia:
            q = gmul(q, gconj(g[A]))
    return q

def Gsep(g, b, a, c):
    """prod n_A over cuts separating b from {a,c}"""
    G = 1
    for A in CUTS:
        ib, ia, ic = (b in A), (a in A), (c in A)
        if ia == ic and ib != ia:
            G *= gnorm(g[A])
    return G

def check_one(rng, r, lo, hi, log):
    primes = random_split_primes(r, lo, hi, rng)
    eps, prof, zeta = build_tuple(primes, rng)
    ok = {}
    gz = ggcd_list(zeta)
    ok['primitive'] = gnorm(gz) == 1
    ok['distinct'] = len(set(zeta)) == 5
    g, chains = gammas_from_profile(primes, prof)
    ok['chain'] = all(chain_ok(chains[p]) for p in primes)
    ok['cox_rep'] = zeta_from_gammas(eps, g) == zeta
    N = gnorm(zeta[0])
    ok['same_norm'] = all(gnorm(z) == N for z in zeta)
    ok['N_is_prod_nA'] = N == eval('*'.join(str(gnorm(g[A])) for A in CUTS))
    # Pi_a
    Pi = {a: gprod([g[A] for A in CUTS if a in A]) for a in [1, 2, 3, 4]}
    ok['ratio'] = all(gmul(zeta[a-1], gconj(Pi[a])) == gmul(zeta[4], Pi[a]) for a in [1, 2, 3, 4])
    ok['dist_a5'] = all(gnorm(gsub(zeta[a-1], zeta[4])) * gnorm(Pi[a]) == 4 * N * Pi[a][1]**2 for a in [1, 2, 3, 4])
    # pairs
    pr = True; cab = True
    for a, b in combinations(FIVE, 2):
        Q = Qpair(g, a, b)
        dab = gnorm(Q)  # = prod n_A over separating cuts
        t = Q[1]
        if gnorm(gsub(zeta[a-1], zeta[b-1])) * dab != 4 * N * t * t:
            pr = False
        gg = ggcd(zeta[a-1], zeta[b-1])
        c = gexactdiv(gsub(zeta[a-1], zeta[b-1]), gg)
        if gnorm(c) != 4 * t * t:
            cab = False
        if a < 5 and b == 5 and Q != Pi[a]:
            pr = False
    ok['pair_t_is_ImQ'] = pr
    ok['c_ab_norm_4t2'] = cab
    # triangle relations over all ordered triples of [5]
    tri = True
    for a, b, c in permutations(FIVE, 3):
        lhs = gmul(Qpair(g, a, b), Qpair(g, b, c))
        G = Gsep(g, b, a, c)
        rhs = Qpair(g, a, c)
        if lhs != (G * rhs[0], G * rhs[1]):
            tri = False
    ok['triangle'] = tri
    # Im/Re/quadratic-form identities for a,b in [4]
    qf = True
    for a, b in combinations([1, 2, 3, 4], 2):
        xa, ta = Pi[a]; xb, tb = Pi[b]
        Gab = 1; Ma = 1; Mb = 1
        for A in CUTS:
            if a in A and b in A: Gab *= gnorm(g[A])
            elif a in A: Ma *= gnorm(g[A])
            elif b in A: Mb *= gnorm(g[A])
        Q = Qpair(g, a, b); y, tab = Q
        if xb*ta - xa*tb != Gab*tab: qf = False
        if xa*xb + ta*tb != Gab*y: qf = False
        if Mb*ta*ta - 2*y*ta*tb + Ma*tb*tb != Gab*tab*tab: qf = False
    ok['Im_Re_quadform'] = qf
    allok = all(ok.values())
    if not allok:
        log.write('FAIL %s %s\n' % (primes, ok))
    return allok, ok

def converse_tests(rng, trials, log):
    """random gamma_A from split primes: (Ch) satisfied <=> primitive."""
    agree = 0; ch_true = 0; tot = 0
    primes = random_split_primes(6, 5, 200, rng)
    for _ in range(trials):
        g = {A: (1, 0) for A in CUTS}
        chains = {p: ([], []) for p in primes}
        for p in primes:
            pi = gauss_prime(p)
            for _k in range(rng.randint(0, 3)):
                A = rng.choice(CUTS)
                if rng.random() < 0.5:
                    g[A] = gmul(g[A], pi); chains[p][0].append(A)
                else:
                    g[A] = gmul(g[A], gconj(pi)); chains[p][1].append(A)
        zs = zeta_from_gammas((1, 0), g)
        prim = gnorm(ggcd_list(zs)) == 1
        ch = all(chain_ok(chains[p]) for p in primes)
        canon = False
        if prim:
            # recover the profile and the canonical gammas
            prof = {}
            for p in primes:
                pi = gauss_prime(p)
                a = [uniquely_divides_count(z, pi) for z in zs]
                e = a[0] + uniquely_divides_count(zs[0], gconj(pi))
                prof[p] = (e, a)
            g2, _ = gammas_from_profile(primes, prof)
            canon = (g2 == g)
        tot += 1; ch_true += ch
        agree += (ch == (prim and canon))
        if ch and not prim:
            log.write('ERROR: (Ch) but imprimitive\n')
        prim_noncanon = locals().get('prim_noncanon', 0)
    log.write('converse: %d/%d trials agree [(Ch) <=> (primitive and gammas canonical)]; (Ch) held in %d\n' % (agree, tot, ch_true))
    return agree == tot

def unit_class_test(rng, trials, log):
    """mix units: multiply zeta_j by i^{k_j}; count pairs where |c_ab|^2/4 is not a perfect square."""
    bad = 0; tot = 0
    for _ in range(trials):
        primes = random_split_primes(5, 5, 300, rng)
        eps, prof, zeta = build_tuple(primes, rng)
        units = [(1, 0), (0, 1), (-1, 0), (0, -1)]
        zz = [gmul(z, rng.choice(units)) for z in zeta]
        for a, b in combinations(range(5), 2):
            if zz[a] == zz[b]:
                continue
            c = gexactdiv(gsub(zz[a], zz[b]), ggcd(zz[a], zz[b]))
            tot += 1
            if gnorm(c) % 4 != 0:
                bad += 1
    log.write('unit-class mixing: %d/%d pairs have |c_ab|^2 not = 0 mod 4 (t_ab not integral)\n' % (bad, tot))

if __name__ == '__main__':
    rng = random.Random(27)
    log = sys.stdout
    n_ok = 0; T = 300
    for k in range(T):
        r = rng.randint(3, 9)
        allok, ok = check_one(rng, r, 5, 400 if k % 2 else 5000, log)
        n_ok += allok
    log.write('cox form: %d/%d random genuine profiles pass all checks %s\n' % (n_ok, T, sorted(ok.keys())))
    converse_tests(rng, 3000, log)
    unit_class_test(rng, 200, log)
