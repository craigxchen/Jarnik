#!/usr/bin/env python3
"""Independent exact check of stab.md Prop 1.1, 1.2 (Gaussian Cox form), Lemma 3.1 (cofactor) and the
chain-deficit example of Prop 3.5.  Own Gaussian arithmetic; Gaussian prime normalisation chosen
DIFFERENT from the author's (here pi_p = a+bi with b > a > 0, author uses a > b > 0)."""
import itertools, random, math, sys
from fractions import Fraction

# ---------- Gaussian integers as (re, im) ----------
def gm(x, y): return (x[0]*y[0] - x[1]*y[1], x[0]*y[1] + x[1]*y[0])
def gc(x): return (x[0], -x[1])
def gn(x): return x[0]*x[0] + x[1]*x[1]
def gsub(x, y): return (x[0]-y[0], x[1]-y[1])
def gadd(x, y): return (x[0]+y[0], x[1]+y[1])
def gpow(x, k):
    r = (1, 0)
    for _ in range(k): r = gm(r, x)
    return r
def gprod(L):
    r = (1, 0)
    for x in L: r = gm(r, x)
    return r
def gdivexact(x, y):
    """x / y if exact in Z[i], else None"""
    n = gn(y); num = gm(x, gc(y))
    if num[0] % n or num[1] % n: return None
    return (num[0]//n, num[1]//n)
def gval(x, pi):
    v = 0
    while True:
        q = gdivexact(x, pi)
        if q is None: return v
        x = q; v += 1
def ggcd(x, y):
    while y != (0, 0):
        n = gn(y); num = gm(x, gc(y))
        q = (round(Fraction(num[0], n)), round(Fraction(num[1], n)))
        x, y = y, gsub(x, gm(q, y))
    return x

def primes1mod4(lo, hi):
    out = []
    for p in range(max(lo, 5), hi):
        if p % 4 == 1 and all(p % q for q in range(2, int(p**0.5)+1)):
            out.append(p)
    return out
def gauss_prime(p):  # b > a > 0 convention (differs from author)
    for a in range(1, int(p**0.5)+1):
        b2 = p - a*a; b = int(round(b2**0.5))
        if b*b == b2 and b > a: return (a, b)
    raise ValueError

UNITS = [(1, 0), (0, 1), (-1, 0), (0, -1)]
CUTS = [frozenset(A) for r in range(1, 5) for A in itertools.combinations(range(1, 5), r)]  # A subset {1..4}

def tuple_from_profile(prof, eps):
    """prof: dict p -> list a_j(p), j=1..5 ; zeta_j = eps prod pi^{a_j} pibar^{e-a_j}"""
    Z = []
    for j in range(5):
        z = eps
        for p, (a, e, pi) in prof.items():
            z = gm(z, gm(gpow(pi, a[j]), gpow(gc(pi), e - a[j])))
        Z.append(z)
    return Z

def layer_gammas(prof):
    """stab.md §1.1 construction: layer (p,t), S={j: a_j>=t}; 5 not in S -> pi into gamma_S; 5 in S -> pibar into
    gamma_{[5]\\S}.  indices j=1..5 ; a list is 0-based (index 4 = point 5)."""
    gam = {A: (1, 0) for A in CUTS}
    for p, (a, e, pi) in prof.items():
        for t in range(1, e+1):
            S = frozenset(j+1 for j in range(5) if a[j] >= t)
            if 5 not in S:
                gam[S] = gm(gam[S], pi)
            else:
                A = frozenset(range(1, 6)) - S
                gam[A] = gm(gam[A], gc(pi))
    return gam

def zeta_from_gammas(gam, eps):
    Z = []
    for a in range(1, 5):
        z = eps
        for A in CUTS:
            z = gm(z, gam[A] if a in A else gc(gam[A]))
        Z.append(z)
    z5 = eps
    for A in CUTS: z5 = gm(z5, gc(gam[A]))
    Z.append(z5)
    return Z

def Qab(gam, a, b):
    """Q_ab = prod_{A ni a, A not ni b} gamma_A * prod_{A ni b, A not ni a} conj(gamma_A); 5 lies in no A"""
    q = (1, 0)
    for A in CUTS:
        ina, inb = (a in A), (b in A)
        if ina and not inb: q = gm(q, gam[A])
        if inb and not ina: q = gm(q, gc(gam[A]))
    return q

def random_profile(rnd, plist, nprimes, emax):
    prof = {}
    for p in rnd.sample(plist, nprimes):
        e = rnd.randint(1, emax)
        while True:
            a = [rnd.randint(0, e) for _ in range(5)]
            if min(a) == 0 and max(a) == e: break
        prof[p] = (a, e, gauss_prime(p))
    return prof

def main():
    rnd = random.Random(20261029)
    plist = primes1mod4(5, 400)
    nfail = {k: 0 for k in ['param', 'primitive', 'norms', 'Qprim', 'zeta_eq_GQ', 'dist', 'tri', 'tri5', 'x_t', 'quad', 'coincide']}
    ncases = 0
    for trial in range(400):
        prof = random_profile(rnd, plist, rnd.randint(2, 7), rnd.choice([1, 1, 2, 3]))
        eps = rnd.choice(UNITS)
        Z = tuple_from_profile(prof, eps)
        g = (0, 0)
        for z in Z: g = ggcd(g, z)
        if gn(g) != 1: nfail['primitive'] += 1
        gam = layer_gammas(prof)
        Z2 = zeta_from_gammas(gam, eps)
        if Z2 != Z: nfail['param'] += 1
        Nprime = gn(Z[0])
        n = {A: gn(gam[A]) for A in CUTS}
        if math.prod(n.values()) != Nprime: nfail['norms'] += 1
        ncases += 1
        Q = {}
        for a in range(1, 6):
            for b in range(1, 6):
                if a == b: continue
                Q[(a, b)] = Qab(gam, a, b)
        for (a, b), q in Q.items():
            # conjugate-primitive: gcd(Q, conj Q) = 1
            if gn(ggcd(q, gc(q))) != 1: nfail['Qprim'] += 1
            # zeta_a = G' Q_ab, zeta_b = G' conj(Q_ab)
            Gp = eps
            for A in CUTS:
                if a in A and b in A: Gp = gm(Gp, gam[A])
                if a not in A and b not in A: Gp = gm(Gp, gc(gam[A]))
            if gm(Gp, q) != Z[a-1] or gm(Gp, gc(q)) != Z[b-1]: nfail['zeta_eq_GQ'] += 1
            # |zeta_a - zeta_b|^2 * N(Q_ab) = 4 t^2 N'
            t = q[1]
            if gn(gsub(Z[a-1], Z[b-1])) * gn(q) != 4 * t * t * Nprime: nfail['dist'] += 1
            if t == 0 and a < b: nfail['coincide'] += 1
        # triangle relation Q_ab Q_bc = G^{(b)}_{abc} Q_ac
        for a, b, c in itertools.permutations(range(1, 6), 3):
            G = 1
            for A in CUTS:
                side = lambda j: (j in A)
                if side(b) != side(a) and side(a) == side(c): G *= n[A]
            if gm(Q[(a, b)], Q[(b, c)]) != (G * Q[(a, c)][0], G * Q[(a, c)][1]): nfail['tri'] += 1
        # Pi_a Pibar_b = G_ab Q_ab ; x_b t_a - x_a t_b = G_ab t_ab ; quadratic identity
        for a, b in itertools.permutations(range(1, 5), 2):
            Pa, Pb = Q[(a, 5)], Q[(b, 5)]
            Gab = math.prod(n[A] for A in CUTS if a in A and b in A)
            if gm(Pa, gc(Pb)) != (Gab * Q[(a, b)][0], Gab * Q[(a, b)][1]): nfail['tri5'] += 1
            xa, ta = Pa; xb, tb = Pb; y, tab = Q[(a, b)]
            if xb * ta - xa * tb != Gab * tab: nfail['x_t'] += 1
            Ma = math.prod(n[A] for A in CUTS if a in A and b not in A)
            Mb = math.prod(n[A] for A in CUTS if b in A and a not in A)
            if Mb * ta * ta - 2 * y * ta * tb + Ma * tb * tb != Gab * tab * tab: nfail['quad'] += 1
    print('Prop 1.1/1.2 on', ncases, 'random genuine one-class profiles (pi=a+bi, b>a>0):')
    for k, v in nfail.items(): print(f'   {k:12s} failures/occurrences: {v}')

    # ---------- (Ch) <=> primitive & canonical on random gamma-vectors ----------
    agree = disagree = 0; ex = None; n_ch = 0
    small = plist[:4]
    for trial in range(3000):
        gam = {A: (1, 0) for A in CUTS}
        for p in rnd.sample(small, rnd.randint(1, 3)):
            pi = gauss_prime(p)
            for _ in range(rnd.randint(1, 4)):
                A = rnd.choice(CUTS)
                gam[A] = gm(gam[A], pi if rnd.random() < 0.5 else gc(pi))
        # chain condition
        ok = True
        for p in small:
            pi = gauss_prime(p)
            for f in (pi, gc(pi)):
                pass
            Spi = [A for A in CUTS for _ in range(gval(gam[A], pi))]
            Spb = [A for A in CUTS for _ in range(gval(gam[A], gc(pi)))]
            for L in (Spi, Spb):
                for X, Y in itertools.combinations(L, 2):
                    if not (X <= Y or Y <= X): ok = False
            u1 = frozenset().union(*Spi) if Spi else frozenset()
            u2 = frozenset().union(*Spb) if Spb else frozenset()
            if u1 & u2: ok = False
        n_ch += ok
        Z = zeta_from_gammas(gam, (1, 0))
        g = (0, 0)
        for z in Z: g = ggcd(g, z)
        prim = gn(g) == 1
        canon = False
        if prim:
            prof = {}
            for p in small:
                pi = gauss_prime(p)
                a = [gval(z, pi) for z in Z]
                e = gval(gn(Z[0]) and Z[0], pi) + gval(Z[0], gc(pi))
                if max(a) > 0 or e > 0:
                    prof[p] = (a, e, pi)
            canon = (layer_gammas(prof) == gam)
        if ok == (prim and canon): agree += 1
        else:
            disagree += 1
            ex = ex or (gam, prim, canon, ok)
    print(f'(Ch) <=> primitive & canonical: agree {agree}, disagree {disagree}; (Ch) held in {n_ch}/3000')
    if ex: print('  example disagreement:', ex)

    # ---------- Prop 1.1(iii) counterexample ----------
    p = 5; pi = gauss_prime(5); pb = gc(pi)
    gam = {A: (1, 0) for A in CUTS}
    for A in (frozenset({1, 3}), frozenset({3}), frozenset({3, 4})): gam[A] = gm(gam[A], pb)
    Z = zeta_from_gammas(gam, (1, 0))
    g = (0, 0)
    for z in Z: g = ggcd(g, z)
    a = [gval(z, pi) for z in Z]
    prof = {5: (a, 3, pi)}
    can = layer_gammas(prof)
    print('Prop 1.1(iii) example: primitive', gn(g) == 1, 'exponents of pi', a,
          'canonical nontrivial blocks', {tuple(sorted(A)): v for A, v in can.items() if v != (1, 0)},
          'vs input', {tuple(sorted(A)): v for A, v in gam.items() if v != (1, 0)})

if __name__ == '__main__':
    main()
