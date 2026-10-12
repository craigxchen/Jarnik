#!/usr/bin/env python3
"""R3. Sections 3-5 (and the algebra of 6.1, 8.1, 8.3) of stability.md on actual clusters.

Circles: N = (content) * prod p^e over small split primes, N <= NMAX, with many lattice points
(r2(N) >= R2MIN), plus content multipliers 2, 9, 49, 4*9, 5^2 (split content).  Every angular window of
M >= 4 consecutive points with chord constant C = maxchord / N^(1/4) <= CMAX is tested.

Exact integer tests:
  T3.1   every prime q | N: exists s in [ceil(M/2), M] with q^(s^2) e^(8F*(s)) N^(s(s-2)) <= X^(2s(s-1)),
         X = max squared chord of the window (this is log q <= max_s b_C(s));
  4.1b   per prime: sum_{x<y} |a_x - a_y| = e M^2/4 - sum_tau (S_tau - M/2)^2  (times 4, integers);
  6.1    Pi = prod (u_x/conj u_x)^(c_x) equals +- prod (pi_p/conj pi_p)^(2 v_p) (cross-multiplied);
  8.1    A_w beta_xy + A_x beta_yw = A_y beta_xw, A's pairwise coprime;
  8.3    a_x + a_y = a_w + a_v  ==>  u_x u_y = u_w u_v exactly (unit 1).
mpmath (40 digits, tolerance 1e-25) tests:
  T3.1 crude forms; T3.3 (1)-(4) with beta' = max_s b'_C(s) (statement) and with C' (stronger);
  Prop 3.5 for several Y; Prop 4.1 (a) lower bound, (c), (4.1), (4.2), (e) incl. Markov counts;
  Prop 5.1 identity and 0 <= X_chord, X_coll, X_prof <= Xi <= (M^2/8)(W'/M - 3log M + K');
  Thm 5.2; Lemma 8.1 norm bounds; |Pi - 1| <= ||c||_1 * (angular width).
Usage: python3 r3_window_scan.py NMAX CMAX R2MIN MMAX
"""
import sys
import random
from math import atan2, pi as PI, comb
import mpmath
from rlib import (factorize, lattice_points, gsub, gmul, gconj, gnorm, ggcd, ggcd_list, gdivexact,
                  split_prime, vpi, gpow, sieve, F8int, F8vec, E)

mpmath.mp.dps = 40
TOL = mpmath.mpf(10) ** -25
K0 = mpmath.mpf('13.9092')
SPLIT = [p for p in sieve(400) if p % 4 == 1]
SMALLP = sieve(10 ** 5)


def F8mp(s):
    return sum(n * mpmath.log(q) for q, n in F8vec(s).items())


def gen_norms(NMAX, R2MIN):
    out = []

    def rec(i, n, r2):
        if r2 >= R2MIN:
            out.append(n)
        for j in range(i, len(SPLIT)):
            p = SPLIT[j]
            m, e = n * p, 1
            if m > NMAX:
                break
            while m <= NMAX:
                rec(j + 1, m, r2 * (e + 1))
                m *= p
                e += 1
    rec(0, 1, 4)
    return sorted(out)


def factor_small(n):
    """factor |beta|^2-size numbers (trial division by SMALLP, remainder prime if < 1e10)."""
    f = {}
    for p in SMALLP:
        if p * p > n:
            break
        while n % p == 0:
            f[p] = f.get(p, 0) + 1
            n //= p
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


def vq_beta(beta, q):
    """largest a with q^a | beta (rational q, Gaussian beta)."""
    a = 0
    x, y = beta
    while x % q == 0 and y % q == 0:
        x //= q
        y //= q
        a += 1
    return a


def chi4(q):
    return 1 if q % 4 == 1 else -1


class Stats(dict):
    def bad(self, key, info=None):
        self[key] = self.get(key, 0) + 1
        if self[key] <= 3:
            print("   FAIL", key, info)


def test_window(Z, N, fac, st, rng):
    M = len(Z)
    W = mpmath.log(N)
    X = max(gnorm(gsub(Z[a], Z[b])) for a in range(M) for b in range(a))
    C = mpmath.sqrt(X) / mpmath.mpf(N) ** mpmath.mpf('0.25')
    logC = mpmath.log(C)
    lCp = max(logC, 0)
    Kp = K0 + 4 * lCp
    K = comb(M, 2)
    st['windows'] = st.get('windows', 0) + 1
    # ---------------- Theorem 3.1 (exact, refined) and crude forms
    F8 = {s: F8int(s) for s in range((M + 1) // 2, M + 1)}
    for q in fac:
        okq = any(q ** (s * s) * F8[s] * N ** (s * (s - 2)) <= X ** (2 * s * (s - 1)) for s in F8)
        if not okq:
            st.bad('T3.1_refined', (N, M, q))
        lq = mpmath.log(q)
        if not lq < 2 * W / M + 4 * lCp:
            st.bad('T3.1_crude1', (N, M, q))
        if not lq <= 2 * W / M - 3 * mpmath.log(mpmath.mpf(M) / 2) + Kp:
            st.bad('T3.1_crude2', (N, M, q))
        r = lq / max(W / s - F8mp(s) / s ** 2 + 4 * (s - 1) / mpmath.mpf(s) * logC for s in F8)
        if r > st.get('T3.1_maxratio', (0,))[0]:
            st['T3.1_maxratio'] = (float(r), N, M, float(C), q)
    # ---------------- normalise
    G = ggcd_list(Z)
    nG = gnorm(G)
    U = [gdivexact(z, G) for z in Z]
    Np = N // nG
    facp = factorize(Np)
    Wp = mpmath.log(Np)
    logCprime = logC - mpmath.log(nG) / 4
    primes = sorted(facp)
    lev = {p: [vpi(u, split_prime(p)) for u in U] for p in primes}
    cuts = [(p, tau) for p in primes for tau in range(1, facp[p] + 1)]
    r = len(cuts)

    def bprime(s, lc):
        return Wp / s + 4 * mpmath.mpf(s - 1) / s * lc - F8mp(s) / s ** 2
    beta_stmt = max(bprime(s, logC) for s in range((M + 1) // 2, M + 1))
    beta_strong = max(bprime(s, logCprime) for s in range((M + 1) // 2, M + 1))
    # ---------------- Theorem 3.3
    classes = {}
    for p, tau in cuts:
        side = frozenset(x for x in range(M) if lev[p][x] < tau)
        key = min(side, frozenset(range(M)) - side, key=lambda S: sorted(S))
        classes[key] = classes.get(key, 0) + mpmath.log(p)
    for wt in classes.values():
        if wt > beta_strong + TOL:
            st.bad('T3.3(1)', (N, M))
    if len(classes) < Wp / beta_stmt - TOL:
        st.bad('T3.3(1)count', (N, M))
    for p in primes:
        if facp[p] * mpmath.log(p) > 2 * beta_strong + TOL:
            st.bad('T3.3(2)', (N, M, p))
        if not facp[p] * mpmath.log(p) < 4 * Wp / M + 8 * lCp:
            st.bad('T3.3(2)crude', (N, M, p))
    om = len(primes)
    if om < Wp / (2 * beta_stmt) - TOL or not om > Wp / (4 * Wp / M + 8 * lCp):
        st.bad('T3.3(3)', (N, M))
    if C <= 1 and not (om > M / 4.0 and len(classes) > M / 2.0 and r > M / 2.0):
        st.bad('T3.3(3)C<=1', (N, M))
    for x in range(M):
        for y in range(x):
            nsep = sum(1 for p, tau in cuts if (lev[p][x] < tau) != (lev[p][y] < tau))
            npr = sum(1 for p in primes if lev[p][x] != lev[p][y])
            bnd = (Wp / 2 + mpmath.log(2) - 2 * lCp)
            if not nsep > bnd / beta_stmt - TOL or not npr > bnd / (2 * beta_stmt) - TOL:
                st.bad('T3.3(4)', (N, M, x, y))
            if C <= 1 and not (nsep > M / 4.0 and npr > M / 8.0):
                st.bad('T3.3(4)C<=1', (N, M))
    # ---------------- Prop 3.5
    for Y in (5, 13, 30, 100, M):
        piY = sum(1 for p in SMALLP if p <= Y)
        WY = sum(facp[p] * mpmath.log(p) for p in primes if p <= Y)
        if not WY < piY * (4 * Wp / M + 8 * lCp):
            st.bad('P3.5a', (N, M, Y))
        if Y >= 2 and not r < Wp / mpmath.log(Y) + piY * (4 * Wp / M + 8 * lCp) / mpmath.log(5):
            st.bad('P3.5b', (N, M, Y))
    # ---------------- pair data, Prop 4.1
    Q4 = 0     # 4Q / log p summed: Q = sum_p log p * sum_tau (S - M/2)^2
    Q = mpmath.mpf(0)
    for p in primes:
        e = facp[p]
        sq4 = sum((2 * sum(1 for a in lev[p] if a < tau) - M) ** 2 for tau in range(1, e + 1))
        dsum = sum(abs(lev[p][x] - lev[p][y]) for x in range(M) for y in range(x))
        if 4 * dsum != e * M * M - sq4:
            st.bad('4.1b_cut_identity', (N, M, p))
        Q += mpmath.log(p) * mpmath.mpf(sq4) / 4
    beta = {}
    gxy = {}
    margins = []
    lb2sum = mpmath.mpf(0)
    for x in range(M):
        for y in range(x):
            d = gsub(U[x], U[y])
            g = ggcd(U[x], U[y])
            b = gdivexact(d, g)
            beta[(x, y)] = b
            gxy[(x, y)] = g
            dxy = sum(abs(lev[p][x] - lev[p][y]) * mpmath.log(p) for p in primes)
            m = dxy - Wp / 2
            margins.append(m)
            Cxy = mpmath.sqrt(gnorm(d)) / mpmath.mpf(Np) ** mpmath.mpf('0.25')
            lb = mpmath.log(gnorm(b)) / 2
            lb2sum += 2 * lb
            if abs(m - (2 * lb - 2 * mpmath.log(Cxy))) > TOL:
                st.bad('4.1a_identity', (N, M))
            odd = sum(vq_beta(b, q) * mpmath.log(q) for q in factor_small(gnorm(b)) if q > 2)
            if m < mpmath.log(2) + 2 * odd - 2 * logCprime - TOL:
                st.bad('4.1a_lower', (N, M))
    msum = sum(margins)
    if abs(msum - (M * Wp / 4 - Q)) > TOL * 10 ** 3:
        st.bad('4.1b_sum', (N, M))
    ls = [l for l in SMALLP if l % 4 == 3 and l + 1 < M]
    Eterm = sum(2 * E(M, l + 1) * mpmath.log(l) for l in ls)
    if Q > M * Wp / 4 + 2 * K * logCprime - K * mpmath.log(2) - Eterm + TOL:
        st.bad('4.1c', (N, M))
    if not Q < M * Wp / 4 - mpmath.mpf(M * M) / 2 * mpmath.log(M) + mpmath.mpf('3.9626') * M * M + 2 * K * lCp:
        st.bad('(4.1)', (N, M))
    mbar = msum / K
    lo = mpmath.mpf(M) / (M - 1) * (mpmath.log(M) - mpmath.mpf('7.9252')) + mpmath.log(2) - 2 * lCp
    if not (lo < mbar <= Wp / (2 * (M - 1)) + TOL):
        st.bad('(4.2)', (N, M))
    if lb2sum > M * Wp / 4 + 2 * K * logCprime - Q + TOL:
        st.bad('4.1e_upper', (N, M))
    if not lb2sum > mpmath.mpf(M * M) / 2 * mpmath.log(M) - mpmath.mpf('3.9626') * M * M + K * mpmath.log(2):
        st.bad('4.1e_lower', (N, M))
    for lam in (0.5, 1, 2):
        big = sum(1 for b in beta.values() if mpmath.log(gnorm(b)) > lam * mpmath.log(M))
        if big > K * (Wp / (2 * (M - 1)) + 2 * lCp) / (lam * mpmath.log(M)) + TOL:
            st.bad('4.1e_markov', (N, M, lam))
    # ---------------- Prop 5.1 decomposition
    Xchord = mpmath.mpf(0)
    kap = {}      # (q, a) -> count
    for (x, y), b in beta.items():
        fb = factor_small(gnorm(b))
        tot = mpmath.log(2) / 2
        for q in fb:
            if q == 2:
                continue
            v = vq_beta(b, q)
            tot += v * mpmath.log(q)
            for a in range(1, v + 1):
                kap[(q, a)] = kap.get((q, a), 0) + 1
        Xchord += mpmath.log(gnorm(b)) / 2 - tot
    # E_{q^a}: for q not dividing N' use E(M, fibre); for p | N' level-wise
    Esum = mpmath.mpf(0)
    qs = set(q for q, a in kap) | set(q for q in SMALLP if q <= 2 * M and q > 2)
    for q in qs:
        a = 1
        while True:
            if q in facp:
                ns = [sum(1 for t in lev[q] if t == lv) for lv in range(facp[q] + 1)]
                Eq = sum(E(n, (q - 1) * q ** (a - 1)) for n in ns)
                fib = (q - 1) * q ** (a - 1)
            else:
                fib = (q - chi4(q)) * q ** (a - 1)
                Eq = E(M, fib)
            if kap.get((q, a), 0) < Eq:
                st.bad('5.1_pigeonhole', (N, M, q, a))
            Esum += Eq * mpmath.log(q)
            if fib >= M and kap.get((q, a), 0) == 0:
                break
            a += 1
    kapsum = sum(c * mpmath.log(q) for (q, a), c in kap.items())
    Xcoll = kapsum - Esum
    lhs51 = K * mpmath.log(2) / 2 + Esum + Q / 2
    Fst = F8mp(M) / 8
    Xprof = lhs51 - Fst
    total = lb2sum / 2 + Q / 2
    if abs(total - (Fst + Xchord + Xcoll + Xprof)) > TOL * 10 ** 3:
        st.bad('5.1_identity', (N, M))
    Xi = M * Wp / 8 + K * logCprime - Fst
    for nm, val in (('Xchord', Xchord), ('Xcoll', Xcoll), ('Xprof', Xprof)):
        if val < -TOL or val > Xi + TOL:
            st.bad('5.1_' + nm, (N, M, float(val), float(Xi)))
    if Xi > mpmath.mpf(M * M) / 8 * (Wp / M - 3 * mpmath.log(M) + Kp) + TOL:
        st.bad('(5.3)', (N, M))
    # ---------------- Theorem 5.2
    lhs52 = sum(mpmath.log(p) / (p - 1) for p in SMALLP if p <= M and p % 4 == 1 and p not in facp)
    if lhs52 > (Wp - 3 * M * mpmath.log(M) + Kp * M) / (2 * M) + TOL:
        st.bad('T5.2', (N, M))
    # ---------------- Prop 6.1 algebra (random zero-sum c) and |Pi - 1| bound
    ang = [mpmath.atan2(z[1], z[0]) for z in U]
    ref = ang[0]
    rel = [mpmath.mpf((a - ref + mpmath.pi) % (2 * mpmath.pi)) - mpmath.pi for a in ang]
    width = max(rel) - min(rel)
    for _ in range(3):
        c = [rng.randint(-2, 2) for _ in range(M)]
        c[-1] -= sum(c)
        if all(t == 0 for t in c):
            continue
        num, den = (1, 0), (1, 0)
        for x in range(M):
            if c[x] > 0:
                num = gmul(num, gpow(U[x], c[x]))
                den = gmul(den, gpow(gconj(U[x]), c[x]))
            elif c[x] < 0:
                num = gmul(num, gpow(gconj(U[x]), -c[x]))
                den = gmul(den, gpow(U[x], -c[x]))
        numR, denR = (1, 0), (1, 0)
        v = {p: sum(c[x] * lev[p][x] for x in range(M)) for p in primes}
        for p in primes:
            pi_ = split_prime(p)
            if v[p] > 0:
                numR = gmul(numR, gpow(pi_, 2 * v[p]))
                denR = gmul(denR, gpow(gconj(pi_), 2 * v[p]))
            elif v[p] < 0:
                numR = gmul(numR, gpow(gconj(pi_), -2 * v[p]))
                denR = gmul(denR, gpow(pi_, -2 * v[p]))
        A1 = gmul(num, denR)
        B1 = gmul(den, numR)
        if not (A1 == B1 or A1 == (-B1[0], -B1[1])):
            st.bad('6.1_Pi_identity', (N, M))
        if all(v[p] == 0 for p in primes) and A1 != B1 and False:
            pass
        Pi = mpmath.mpc(*num) / mpmath.mpc(*den)
        if abs(Pi - 1) > sum(abs(t) for t in c) * width + TOL:
            st.bad('6.1_Pi_bound', (N, M))
        if any(v[p] != 0 for p in primes) and abs(Pi - 1) < TOL:
            st.bad('6.1_Pi_ne_1', (N, M))
        st['6.1_tests'] = st.get('6.1_tests', 0) + 1
    # ---------------- Lemma 8.1 triangles, Prop 8.3 parallel chords
    for _ in range(min(6, comb(M, 3))):
        x, y, w = sorted(rng.sample(range(M), 3), reverse=True)
        G3 = ggcd_list([U[x], U[y], U[w]])
        Aw = gdivexact(gxy[(x, y)], G3)
        Ax = gdivexact(gxy[(y, w)], G3)
        Ay = gdivexact(gxy[(x, w)], G3)
        lhs = gsub(gmul(Ay, gmul(gxy[(x, w)], beta[(x, w)])), (0, 0))
        # identity: A_w beta_xy + A_x beta_yw = A_y beta_xw
        t1 = gmul(Aw, beta[(x, y)])
        t2 = gmul(Ax, beta[(y, w)])
        t3 = gmul(Ay, beta[(x, w)])
        if (t1[0] + t2[0], t1[1] + t2[1]) != t3:
            st.bad('8.1_identity', (N, M))
        for P1, P2 in ((Aw, Ax), (Aw, Ay), (Ax, Ay)):
            if gnorm(ggcd(P1, P2)) != 1:
                st.bad('8.1_coprime', (N, M))
        g3 = mpmath.log(gnorm(G3))
        mxy = margins[[(a, b) for a in range(M) for b in range(a)].index((x, y))]
        la = mpmath.log(gnorm(Aw))
        if abs(la - (Wp / 2 - mxy - g3)) > TOL * 100:
            st.bad('8.1_norm', (N, M))
        if la < Wp / 6 - mxy - mpmath.mpf(8) / 3 * lCp + mpmath.mpf(4) / 3 * mpmath.log(2) - TOL:
            st.bad('8.1_lower', (N, M))
        if la > Wp / 2 - mpmath.log(2) + 2 * lCp + TOL:
            st.bad('8.1_upper', (N, M))
        st['8.1_tests'] = st.get('8.1_tests', 0) + 1
    for _ in range(10):
        if M < 4:
            break
        x, y, w, v = rng.sample(range(M), 4)
        if all(lev[p][x] + lev[p][y] == lev[p][w] + lev[p][v] for p in primes):
            st['8.3_parallel'] = st.get('8.3_parallel', 0) + 1
            if gmul(U[x], U[y]) != gmul(U[w], U[v]):
                st.bad('8.3_unit', (N, M))


def main():
    NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 10 ** 10
    CMAX = float(sys.argv[2]) if len(sys.argv) > 2 else 8.0
    R2MIN = int(sys.argv[3]) if len(sys.argv) > 3 else 256
    MMAX = int(sys.argv[4]) if len(sys.argv) > 4 else 14
    rng = random.Random(2026)
    base = gen_norms(NMAX, R2MIN)
    norms = []
    for n in base:
        norms.append(n)
    extra = []
    for n in base[::7]:
        for c in (2, 9, 49, 36, 25):
            if n * c <= 4 * NMAX:
                extra.append(n * c)
    norms += extra
    print(f"circles: {len(base)} base norms <= {NMAX} with r2 >= {R2MIN}, plus {len(extra)} with content")
    st = Stats()
    Mhist = {}
    for N in norms:
        fac = factorize(N)
        pts = sorted(lattice_points(fac), key=lambda z: atan2(z[1], z[0]) % (2 * PI))
        n = len(pts)
        q4 = N ** 0.25
        for i in range(n):
            for M in range(4, min(MMAX, n // 2) + 1):
                Z = [pts[(i + t) % n] for t in range(M)]
                X = max(gnorm(gsub(Z[a], Z[b])) for a in range(M) for b in range(a))
                if X ** 0.5 / q4 > CMAX:
                    break
                Mhist[M] = Mhist.get(M, 0) + 1
                test_window(Z, N, fac, st, rng)
    print("windows per M:", dict(sorted(Mhist.items())))
    fails = {k: v for k, v in st.items() if not k.endswith('tests') and k not in
             ('windows', 'T3.1_maxratio', '8.3_parallel')}
    print("stats:", {k: v for k, v in st.items() if k.endswith('tests') or k in ('windows', '8.3_parallel')})
    print("tightest T3.1 ratio log q / max_s b_C(s):", st.get('T3.1_maxratio'))
    print("failures:", fails)
    print("ALL CHECKS PASSED:", len(fails) == 0)


if __name__ == "__main__":
    main()
