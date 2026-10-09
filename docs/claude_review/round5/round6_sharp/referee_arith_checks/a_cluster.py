"""Referee check (arith lens) of sharp.md Prop. 1.2 (actual clusters lie in F^all_X) on genuine
lattice-point clusters, with exact Gaussian-integer arithmetic.  Independent code.

For several N (products of primes = 1 mod 4, some with exponent 2) and several windows of M
consecutive lattice points of x^2 + y^2 = N (in angle), after gcd normalisation:
 K1  (1.1) of outside.md: integer windings k_x with <2a_x - e, phi> = theta* + delta_x + (pi/2)k_x,
     and z_x = eps_0 i^(-k_x) prod pi^a conj(pi)^(e-a) with ONE common unit eps_0;
 K2  r_j(q) = pi_j mod q^(A_q): Norm r_j = p_j, rho_x(q) = conj(eps_0) z_x mod q^(A_q);
 K3  (e-char): for every c with ||c||_1 <= NMAX, sum c = 0, v(c) != 0 and every eta: with
     G = gcd(U, V), D = U/G - eta V/G:  D != 0,  (1+i) | D,  q^(a_(q,eta)(c)) | D for every q,
     Norm(G) * prod p^|v_j| = N^(n/2), U/G = unit*A_v, Norm(U - eta V) = Norm(G) Norm(D);
     hence |sin((c.delta - arg eta)/2)| >= 2^(-1/2) m e^(-w/2); the float value of the left side
     agrees with |U - eta V| / (2 N^(n/4));
 K4  (e-all): for every v with ||v||_1 <= VMAX and every s: D_s in {Y, X-Y, X, X+Y} of A_v is
     nonzero, divisible by m_(v,s), and |D_s| nu_s = |A_v| |sin(<v,phi> - s pi/4)|;
 K5  same-level clause: for p_l <= X dividing N, level residues rho_x(p_l) (product over j != l)
     equal conj(eps_0) r_x (two_thirds' r_x) mod p_l^A, and for pairs on a common level,
     p_l^(a) | c_xy = (z_x - z_y)/gcd(z_x, z_y) whenever rho_x = rho_y mod p_l^a;
 K6  Cor. 2.2(4) on genuine data with the actual windings: the class of rho_x/rho_y in T/T^2 equals
     sum_j (a_xj - a_yj) nu_j(q) + lambda_q (k_x - k_y) mod 2.
"""
import math, itertools, random

rng = random.Random(7)


def is_prime(n):
    if n < 2:
        return False
    d = 2
    while d * d <= n:
        if n % d == 0:
            return False
        d += 1
    return True


def gmul(a, b):
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])


def gconj(a):
    return (a[0], -a[1])


def gnorm(a):
    return a[0] * a[0] + a[1] * a[1]


def gpow(a, e):
    r = (1, 0)
    for _ in range(e):
        r = gmul(r, a)
    return r


def gdivmod(a, b):
    # nearest-integer division in Z[i]
    n = gnorm(b)
    num = gmul(a, gconj(b))
    qx = (2 * num[0] + n) // (2 * n)
    qy = (2 * num[1] + n) // (2 * n)
    q = (qx, qy)
    r = (a[0] - gmul(q, b)[0], a[1] - gmul(q, b)[1])
    return q, r


def gdiv_exact(a, b):
    q, r = gdivmod(a, b)
    assert r == (0, 0), (a, b)
    return q


def gdivides(b, a):
    return gdivmod(a, b)[1] == (0, 0)


def ggcd(a, b):
    while b != (0, 0):
        _, r = gdivmod(a, b)
        a, b = b, r
    return a


UNITS = [(1, 0), (0, 1), (-1, 0), (0, -1)]


def two_squares(p):
    for x in range(1, math.isqrt(p) + 1):
        y = math.isqrt(p - x * x)
        if x * x + y * y == p:
            return (x, y)


def mmod(a, Q):
    return (a[0] % Q, a[1] % Q)


def mmul(a, b, Q):
    return ((a[0] * b[0] - a[1] * b[1]) % Q, (a[0] * b[1] + a[1] * b[0]) % Q)


def mpow(a, e, Q):
    if e < 0:
        n = pow(gnorm(a) % Q, -1, Q)
        a = ((a[0] * n) % Q, (-a[1] * n) % Q)
        e = -e
    r, b = (1, 0), mmod(a, Q)
    while e:
        if e & 1:
            r = mmul(r, b, Q)
        b = mmul(b, b, Q); e >>= 1
    return r


def val(z, pi):
    v = 0
    while gdivides(pi, z):
        z = gdiv_exact(z, pi); v += 1
    return v


def leg(a, q):
    a %= q
    if a == 0:
        return 0
    return 1 if pow(a, (q - 1) // 2, q) == 1 else -1


def run(primes_exp, M, start, X, NMAX, VMAX, stats):
    pis = {p: two_squares(p) for p, _ in primes_exp}
    # all lattice points of norm N
    pts = []
    for exps in itertools.product(*[range(e + 1) for _, e in primes_exp]):
        z = (1, 0)
        for (p, e), a in zip(primes_exp, exps):
            z = gmul(z, gmul(gpow(pis[p], a), gpow(gconj(pis[p]), e - a)))
        for u in UNITS:
            pts.append(gmul(u, z))
    pts.sort(key=lambda z: math.atan2(z[1], z[0]) % (2 * math.pi))
    win = [pts[(start + i) % len(pts)] for i in range(M)]
    # normalise
    g = win[0]
    for z in win[1:]:
        g = ggcd(g, z)
    zs = [gdiv_exact(z, g) for z in win]
    N = gnorm(zs[0])
    assert all(gnorm(z) == N for z in zs)
    # exponent data relative to the primes still dividing N
    plist = [p for p, _ in primes_exp if N % p == 0]
    e = [0] * len(plist)
    for j, p in enumerate(plist):
        nn = N
        while nn % p == 0:
            nn //= p; e[j] += 1
    assert math.prod(p ** ej for p, ej in zip(plist, e)) == N
    r = len(plist)
    a = [[val(z, pis[p]) for p in plist] for z in zs]
    ok = {}
    ok['norm'] = all(min(a[x][j] for x in range(M)) == 0 and max(a[x][j] for x in range(M)) == e[j] for j in range(r))
    core = []
    for x in range(M):
        w = (1, 0)
        for j, p in enumerate(plist):
            w = gmul(w, gmul(gpow(pis[p], a[x][j]), gpow(gconj(pis[p]), e[j] - a[x][j])))
        core.append(w)
    phi = [math.atan2(pis[p][1], pis[p][0]) for p in plist]
    ang = [math.atan2(z[1], z[0]) for z in zs]
    th = ang[0]
    delta = [((t - th) % (2 * math.pi)) for t in ang]
    Delta = max(delta)
    ok['short'] = Delta < math.pi / 2
    k = []
    resid = 0
    for x in range(M):
        lhs = sum((2 * a[x][j] - e[j]) * phi[j] for j in range(r))
        kk = (lhs - th - delta[x]) / (math.pi / 2)
        k.append(round(kk)); resid = max(resid, abs(kk - round(kk)))
    ok['K1_integral_k'] = resid < 1e-7
    eps0s = set()
    for x in range(M):
        # z_x = eps0 * i^(-k_x) * core_x
        ik = UNITS[(-k[x]) % 4]
        t = gmul(ik, core[x])
        e0 = gdiv_exact(zs[x], t)
        eps0s.add(e0)
    ok['K1_common_unit'] = len(eps0s) == 1 and list(eps0s)[0] in UNITS
    eps0 = list(eps0s)[0]
    Rr = math.sqrt(N)
    C = (Delta * Rr) / math.sqrt(Rr)   # arc length / sqrt(R)
    W = math.log(N)
    logp = [math.log(p) for p in plist]
    # moduli
    mods = []
    for q in range(3, X + 1, 2):
        if not is_prime(q):
            continue
        A = 1
        while q ** (A + 1) <= X:
            A += 1
        mods.append((q, A))
    good = {kk: True for kk in ['K2', 'K3', 'K4', 'K5', 'K6']}
    resq = {}
    for (q, A) in mods:
        Q = q ** A
        if N % q == 0:
            continue
        rj = [mmod(pis[p], Q) for p in plist]
        good['K2'] &= all((x * x + y * y) % Q == p % Q for (x, y), p in zip(rj, plist))
        rho = []
        for x in range(M):
            v = mpow((0, 1), -k[x], Q)
            for j in range(r):
                v = mmul(v, mmul(mpow(rj[j], a[x][j], Q), mpow(gconj(rj[j]), e[j] - a[x][j], Q), Q), Q)
            rho.append(v)
            good['K2'] &= v == mmod(gmul(gconj(eps0), zs[x]), Q)
        resq[q] = (A, rho, rj)
        # K6: class in T/T^2 of rho_x/rho_y: Legendre of Norm(w) for w with w/conj(w) = rho_x/rho_y
        # use rho_x/rho_y = t, t in T_Q ; nu(t) = Legendre(Norm(1 + t))-trick is fragile; instead
        # compute nu via the square test t^(m/2) = 1
        chi = 1 if q % 4 == 1 else -1
        m = (q - chi) * q ** (A - 1)
        lam = 1 if q % 8 in (3, 5) else 0
        nu = [1 if leg(p, q) == -1 else 0 for p in plist]
        for x in range(M):
            t = mmul(rho[x], mpow(rho[0], -1, Q), Q)
            issq = mpow(t, m // 2, Q) == (1, 0)
            par = (sum((a[x][j] - a[0][j]) * nu[j] for j in range(r)) + lam * (k[x] - k[0])) % 2
            good['K6'] &= (issq == (par == 0))
    # characters
    cs = []
    for n in range(2, NMAX + 1, 2):
        for supp in itertools.combinations(range(M), min(M, n)):
            pass
    allc = [c for c in itertools.product(range(-3, 4), repeat=M) if sum(c) == 0 and 0 < sum(map(abs, c)) <= NMAX]
    if len(allc) > 4000:
        allc = rng.sample(allc, 4000)
    nchar = 0
    worst_ratio = 0.0
    for c in allc:
        v = [sum(c[x] * a[x][j] for x in range(M)) for j in range(r)]
        if not any(v):
            continue
        n = sum(map(abs, c))
        U, V = (1, 0), (1, 0)
        for x in range(M):
            if c[x] > 0:
                U = gmul(U, gpow(zs[x], c[x]))
            elif c[x] < 0:
                V = gmul(V, gpow(zs[x], -c[x]))
        G = ggcd(U, V)
        Av = (1, 0)
        for j, p in enumerate(plist):
            Av = gmul(Av, gpow(pis[p], v[j]) if v[j] > 0 else gpow(gconj(pis[p]), -v[j]))
        UG = gdiv_exact(U, G); VG = gdiv_exact(V, G)
        good['K3'] &= any(gmul(u, Av) == UG for u in UNITS) and any(gmul(u, gconj(Av)) == VG for u in UNITS)
        good['K3'] &= gnorm(G) * math.prod(p ** abs(vj) for p, vj in zip(plist, v)) == N ** (n // 2)
        w = sum(abs(vj) * lp for vj, lp in zip(v, logp))
        cd = sum(c[x] * delta[x] for x in range(M))
        for s, eta in enumerate(UNITS):
            D = (UG[0] - gmul(eta, VG)[0], UG[1] - gmul(eta, VG)[1])
            good['K3'] &= D != (0, 0) and gdivides((1, 1), D)
            mlog = 0.0
            for q, (A, rho, rj) in resq.items():
                aq = 0
                for lev in range(1, A + 1):
                    qa = q ** lev
                    P = (1, 0)
                    for x in range(M):
                        if c[x]:
                            P = mmul(P, mpow(rho[x], c[x], qa), qa)
                    if P == mmod(eta, qa):
                        aq = lev
                    else:
                        break
                if aq:
                    good['K3'] &= D[0] % q ** aq == 0 and D[1] % q ** aq == 0
                    mlog += aq * math.log(q)
            UV = (U[0] - gmul(eta, V)[0], U[1] - gmul(eta, V)[1])
            good['K3'] &= gnorm(UV) == gnorm(G) * gnorm(D)
            # float: |sin((c.delta - arg eta)/2)| vs |U - eta V| / (2 N^(n/4))
            lhs = abs(math.sin((cd - s * math.pi / 2) / 2))
            if lhs > 1e-6:
                rel = math.exp(0.5 * math.log(gnorm(UV)) - math.log(2) - (n / 4) * W)
                good['K3'] &= abs(rel - lhs) <= 1e-6 * lhs + 1e-12
                rhs_log = -0.5 * math.log(2) + mlog - w / 2
                good['K3'] &= math.log(lhs) >= rhs_log - 1e-9
                worst_ratio = max(worst_ratio, rhs_log - math.log(lhs))
            else:
                # tiny: rely on exact identity Norm(UV) = Norm(G) Norm(D) and Norm(D) >= 2 m^2
                good['K3'] &= gnorm(D) >= 2 * math.exp(2 * mlog) * (1 - 1e-12)
        nchar += 1
    # monomials
    nv = 0
    for vv in itertools.product(range(-2, 3), repeat=r):
        if not any(vv) or sum(map(abs, vv)) > VMAX:
            continue
        Av = (1, 0)
        for j, p in enumerate(plist):
            Av = gmul(Av, gpow(pis[p], vv[j]) if vv[j] > 0 else gpow(gconj(pis[p]), -vv[j]))
        Xa, Ya = Av
        Ds = [Ya, Xa - Ya, Xa, Xa + Ya]
        nus = [1.0, 2 ** -0.5, 1.0, 2 ** -0.5]
        argA = sum(vj * ph for vj, ph in zip(vv, phi))
        for s in range(4):
            good['K4'] &= Ds[s] != 0
            for q, (A, rho, rj) in resq.items():
                aq = 0
                for lev in range(1, A + 1):
                    qa = q ** lev
                    P = (1, 0)
                    for j in range(r):
                        if vv[j]:
                            t = mmul(rj[j], mpow(gconj(rj[j]), -1, qa), qa)
                            P = mmul(P, mpow(t, vv[j], qa), qa)
                    if P == mpow((0, 1), s, qa):
                        aq = lev
                    else:
                        break
                if aq:
                    good['K4'] &= Ds[s] % q ** aq == 0
            lhs = math.sqrt(gnorm(Av)) * abs(math.sin(argA - s * math.pi / 4))
            good['K4'] &= abs(abs(Ds[s]) * nus[s] - lhs) <= 1e-7 * max(1.0, lhs)
        nv += 1
    # same-level clause
    nlev = 0
    for l, p in enumerate(plist):
        if p > X:
            continue
        A = 1
        while p ** (A + 1) <= X:
            A += 1
        Q = p ** A
        piL = pis[p]
        for x in range(M):
            rx = gdiv_exact(zs[x], gmul(gpow(piL, a[x][l]), gpow(gconj(piL), e[l] - a[x][l])))
            v = mpow((0, 1), -k[x], Q)
            for j in range(r):
                if j == l:
                    continue
                rj = mmod(pis[plist[j]], Q)
                v = mmul(v, mmul(mpow(rj, a[x][j], Q), mpow(gconj(rj), e[j] - a[x][j], Q), Q), Q)
            good['K5'] &= v == mmod(gmul(gconj(eps0), rx), Q)
            for y in range(x + 1, M):
                if a[x][l] != a[y][l]:
                    continue
                ry = gdiv_exact(zs[y], gmul(gpow(piL, a[y][l]), gpow(gconj(piL), e[l] - a[y][l])))
                lev = 0
                for t in range(1, A + 1):
                    if mmod(rx, p ** t) == mmod(ry, p ** t):
                        lev = t
                cxy = gdiv_exact((zs[x][0] - zs[y][0], zs[x][1] - zs[y][1]), ggcd(zs[x], zs[y]))
                good['K5'] &= cxy[0] % p ** lev == 0 and cxy[1] % p ** lev == 0
                nlev += 1
    stats.append((N, M, r, e, round(C, 3), k, nchar, nv, nlev, worst_ratio))
    ok.update(good)
    return ok


if __name__ == "__main__":
    configs = [
        ([(5, 1), (13, 1), (17, 1), (29, 1), (37, 1), (41, 1)], 6),
        ([(5, 2), (13, 1), (17, 2), (29, 1), (53, 1)], 6),
        ([(5, 1), (13, 2), (37, 1), (41, 1), (61, 1), (73, 1)], 7),
        ([(5, 3), (17, 1), (29, 1), (89, 1), (97, 1)], 5),
        ([(13, 1), (17, 1), (29, 1), (37, 1), (41, 1), (53, 1), (61, 1)], 6),
    ]
    allok = True
    stats = []
    for prs, M in configs:
        npts = 4 * math.prod(e + 1 for _, e in prs)
        for start in rng.sample(range(npts), 3):
            ok = run(prs, M, start, X=60, NMAX=6, VMAX=4, stats=stats)
            if not all(ok.values()):
                print("FAIL", prs, M, start, ok)
            allok &= all(ok.values())
    for s in stats:
        print(f"N={s[0]} M={s[1]} r={s[2]} e={s[3]} C={s[4]} k={s[5]} chars={s[6]} monomials={s[7]} same-level pairs={s[8]} max(rhs-lhs log)={s[9]:.3f}")
    print("ALL ACTUAL-CLUSTER (PROP 1.2) CHECKS PASSED" if allok else "PROP 1.2 FAILURE")
