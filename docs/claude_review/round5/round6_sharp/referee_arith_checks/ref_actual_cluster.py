"""Referee check of sharp.md Prop. 1.2: actual normalised clusters lie in F^char_X and F^all_X.

For several N (products of primes = 1 mod 4, some squared, some <= X so that level data occur) we
take every window of M consecutive lattice points (by angle) on x^2 + y^2 = N, normalise by the
Gaussian gcd (two_thirds Lemma 1.1), compute the profile a_x w.r.t. FIXED Gaussian primes pi_j,
real lifts phi_j = arg pi_j, offsets delta_x, and windings k_x from (1.1) (checked to be integers).
Residue data r_j(q) = pi_j mod q^(A_q) for every odd prime q <= X (level data at q = p_l | N).
Then, with floating-point angles (independent of the exact proof identity):
  (e-char) for every c in [-2,2]^M, sum c = 0, v(c) != 0, every unit eta:
           |sin((c.delta - arg eta)/2)| >= 2^(-1/2) m_(c,eta) e^(-w(v(c))/2),
           m from the point residues rho_x(q) = i^(-k_x) prod r_j^(a_xj) conj(r_j)^(e_j-a_xj),
           plus the same-level clause for pairs (eta = 1) at the p_l <= X;
  (e-all)  for every nonzero v in [-2,2]^r and s in Z/4:
           |sin(<v,phi> - s pi/4)| >= nu_s m_(v,s) e^(-w(v)/2).
Reports the minimal ratio LHS/RHS (must be >= 1) and how often collisions occurred."""
import itertools, math

X = 300


def is_prime(n):
    if n < 2:
        return False
    i = 2
    while i * i <= n:
        if n % i == 0:
            return False
        i += 1
    return True


ODDQ = [q for q in range(3, X + 1) if is_prime(q)]


def gmul(a, b):
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])


def gconj(a):
    return (a[0], -a[1])


def gdivexact(a, b):
    n = b[0] ** 2 + b[1] ** 2
    t = gmul(a, gconj(b))
    if t[0] % n or t[1] % n:
        return None
    return (t[0] // n, t[1] // n)


def gmod(a, b):
    n = b[0] ** 2 + b[1] ** 2
    t = gmul(a, gconj(b))
    qx = (2 * t[0] + n) // (2 * n)
    qy = (2 * t[1] + n) // (2 * n)
    pr = gmul(b, (qx, qy))
    return (a[0] - pr[0], a[1] - pr[1])


def ggcd(a, b):
    while b != (0, 0):
        a, b = b, gmod(a, b)
    return a


def gpow_mod(a, e, Q):
    r = (1, 0)
    b = (a[0] % Q, a[1] % Q)
    if e < 0:
        n = pow((b[0] ** 2 + b[1] ** 2) % Q, -1, Q)
        b = ((b[0] * n) % Q, (-b[1] * n) % Q)
        e = -e
    while e:
        if e & 1:
            r = gmul(r, b)
            r = (r[0] % Q, r[1] % Q)
        b = gmul(b, b)
        b = (b[0] % Q, b[1] % Q)
        e >>= 1
    return r


def two_squares(p):
    for a in range(1, int(math.isqrt(p)) + 1):
        b2 = p - a * a
        b = math.isqrt(b2)
        if b * b == b2:
            return (a, b)


def run(Nfac, Mlist):
    N = 1
    for p, e in Nfac.items():
        N *= p ** e
    pts = []
    for x in range(-math.isqrt(N), math.isqrt(N) + 1):
        y2 = N - x * x
        y = math.isqrt(y2)
        if y * y == y2:
            pts.append((x, y))
            if y:
                pts.append((x, -y))
    pts.sort(key=lambda z: math.atan2(z[1], z[0]))
    stats = dict(tight_char_m_gt_1=0, tight_mon_m_gt_1=0, tight_char=0, tight_mon=0, clusters=0, chars=0, mons=0, coll_char=0, coll_level=0, coll_mon=0,
                 minratio_char=float("inf"), minratio_mon=float("inf"), kmax=0, k_int_err=0.0)
    for M in Mlist:
        cap = 10 if M <= 4 else 3
        done = 0
        for st in range(0, len(pts), max(1, len(pts) // 25)):
            if done >= cap:
                break
            Z = [pts[(st + t) % len(pts)] for t in range(M)]
            ang = [math.atan2(z[1], z[0]) for z in Z]
            spread = (ang[-1] - ang[0]) % (2 * math.pi)
            if spread > 1.0:      # keep it an arc (any arc works for Prop 1.2; avoid wrap issues)
                continue
            g = Z[0]
            for z in Z[1:]:
                g = ggcd(g, z)
            Zn = [gdivexact(z, g) for z in Z]
            Nn = Zn[0][0] ** 2 + Zn[0][1] ** 2
            # factor Nn over the primes of N
            fac = {}
            t = Nn
            for p in Nfac:
                while t % p == 0:
                    fac[p] = fac.get(p, 0) + 1
                    t //= p
            assert t == 1
            primes = sorted(fac)
            r = len(primes)
            pis = [two_squares(p) for p in primes]
            e = [fac[p] for p in primes]
            # exponents a_xj = v_{pi_j}(z_x); unit eps_x
            Amat = []
            for z in Zn:
                row = []
                u = z
                for pi in pis:
                    v = 0
                    while True:
                        w_ = gdivexact(u, pi)
                        if w_ is None:
                            break
                        u = w_
                        v += 1
                    row.append(v)
                # remaining u must be conj(pi)^(e-a) * unit
                for j, pi in enumerate(pis):
                    for _ in range(e[j] - row[j]):
                        u = gdivexact(u, gconj(pi))
                        assert u is not None
                assert u in [(1, 0), (0, 1), (-1, 0), (0, -1)]
                Amat.append(row)
            for j in range(r):  # normalisation (Lemma 1.1(d))
                col = [Amat[x][j] for x in range(M)]
                assert min(col) == 0 and max(col) == e[j]
            stats["clusters"] += 1
            done += 1
            phi = [math.atan2(pi[1], pi[0]) for pi in pis]
            thn = [math.atan2(z[1], z[0]) for z in Zn]
            delta = [(thn[x] - thn[0]) % (2 * math.pi) for x in range(M)]
            delta = [d if d < math.pi else d - 2 * math.pi for d in delta]
            assert min(delta) >= -1e-12   # windows sorted by angle; theta* = first point
            lin = [sum((2 * Amat[x][j] - e[j]) * phi[j] for j in range(r)) for x in range(M)]
            theta0 = lin[0] - delta[0]
            kk = []
            for x in range(M):
                kf = (lin[x] - theta0 - delta[x]) / (math.pi / 2)
                k = round(kf)
                stats["k_int_err"] = max(stats["k_int_err"], abs(kf - k))
                kk.append(k)
            stats["kmax"] = max(stats["kmax"], max(abs(k) for k in kk))
            logp = [math.log(p) for p in primes]
            # residues
            res = {}
            for q in ODDQ:
                Aq = 1
                while q ** (Aq + 1) <= X:
                    Aq += 1
                Q = q ** Aq
                if Nn % q:
                    rho = []
                    for x in range(M):
                        u = gpow_mod((0, 1), (-kk[x]) % 4, Q)
                        for j in range(r):
                            u = gmul(u, gpow_mod(pis[j], Amat[x][j], Q))
                            u = gmul(u, gpow_mod(gconj(pis[j]), e[j] - Amat[x][j], Q))
                            u = (u[0] % Q, u[1] % Q)
                        assert (u[0] ** 2 + u[1] ** 2 - Nn) % Q == 0      # Norm rho_x = N
                        rho.append(u)
                    tj = [gmul(gpow_mod(pis[j], 1, Q), gpow_mod(gconj(pis[j]), -1, Q)) for j in range(r)]
                    tj = [(t_[0] % Q, t_[1] % Q) for t_ in tj]
                    res[q] = ("free", Q, Aq, rho, tj)
                else:
                    l = primes.index(q)
                    rho = []
                    for x in range(M):
                        u = gpow_mod((0, 1), (-kk[x]) % 4, Q)
                        for j in range(r):
                            if j == l:
                                continue
                            u = gmul(u, gpow_mod(pis[j], Amat[x][j], Q))
                            u = gmul(u, gpow_mod(gconj(pis[j]), e[j] - Amat[x][j], Q))
                            u = (u[0] % Q, u[1] % Q)
                        rho.append(u)
                    res[q] = ("level", Q, Aq, rho, l)
            units = [(1, 0), (0, 1), (-1, 0), (0, -1)]

            def coll_level(prod_fn, Q, Aq, q, target):
                a = 0
                for lev in range(1, Aq + 1):
                    Ql = q ** lev
                    u = prod_fn(Ql)
                    if (u[0] - target[0]) % Ql == 0 and (u[1] - target[1]) % Ql == 0:
                        a = lev
                    else:
                        break
                return a

            # (e-char)
            for c in itertools.product(range(-2, 3), repeat=M):
                if sum(c) != 0 or not any(c):
                    continue
                v = [sum(c[x] * Amat[x][j] for x in range(M)) for j in range(r)]
                if not any(v):
                    continue
                w = sum(abs(v[j]) * logp[j] for j in range(r))
                cd = sum(c[x] * delta[x] for x in range(M))
                ispair = sorted(abs(t_) for t_ in c if t_) == [1, 1]
                for s in range(4):
                    eta = units[s]
                    logm = 0.0
                    for q, data in res.items():
                        if data[0] == "free":
                            _, Q, Aq, rho, _ = data

                            def pf(Ql, rho=rho):
                                u = (1, 0)
                                for x in range(M):
                                    if c[x]:
                                        u = gmul(u, gpow_mod(rho[x], c[x], Ql))
                                        u = (u[0] % Ql, u[1] % Ql)
                                return u
                            a = coll_level(pf, Q, Aq, q, eta)
                            if a:
                                stats["coll_char"] += 1
                            logm += a * math.log(q)
                        elif ispair and s == 0:
                            _, Q, Aq, rho, l = data
                            xs = [x for x in range(M) if c[x]]
                            if Amat[xs[0]][l] == Amat[xs[1]][l]:
                                def pf(Ql, rho=rho, xs=xs):
                                    u = gmul(rho[xs[0]], gpow_mod(rho[xs[1]], -1, Ql))
                                    return (u[0] % Ql, u[1] % Ql)
                                a = coll_level(pf, Q, Aq, q, (1, 0))
                                if a:
                                    stats["coll_level"] += 1
                                logm += a * math.log(q)
                    lhs = abs(math.sin((cd - s * math.pi / 2) / 2))
                    logrhs = -0.5 * math.log(2) + logm - w / 2
                    ratio = math.log(lhs) - logrhs if lhs > 0 else -float("inf")
                    stats["minratio_char"] = min(stats["minratio_char"], ratio)
                    stats["tight_char"] += abs(ratio) < 1e-9
                    stats["tight_char_m_gt_1"] += (abs(ratio) < 1e-9) and logm > 0
                    stats["chars"] += 1
            # (e-all)
            for v in itertools.product(range(-2, 3), repeat=r):
                if not any(v):
                    continue
                w = sum(abs(v[j]) * logp[j] for j in range(r))
                ang_v = sum(v[j] * phi[j] for j in range(r))
                for s in range(4):
                    logm = 0.0
                    for q, data in res.items():
                        if data[0] != "free":
                            continue
                        _, Q, Aq, rho, tj = data

                        def pf(Ql, tj=tj):
                            u = (1, 0)
                            for j in range(r):
                                if v[j]:
                                    u = gmul(u, gpow_mod(tj[j], v[j], Ql))
                                    u = (u[0] % Ql, u[1] % Ql)
                            return u
                        a = coll_level(pf, Q, Aq, q, units[s])
                        if a:
                            stats["coll_mon"] += 1
                        logm += a * math.log(q)
                    nus = 1.0 if s % 2 == 0 else 2 ** -0.5
                    lhs = abs(math.sin(ang_v - s * math.pi / 4))
                    ratio = math.log(lhs) - (math.log(nus) + logm - w / 2)
                    stats["minratio_mon"] = min(stats["minratio_mon"], ratio)
                    stats["tight_mon"] += abs(ratio) < 1e-9
                    stats["tight_mon_m_gt_1"] += (abs(ratio) < 1e-9) and logm > 0
                    stats["mons"] += 1
    return stats


allok = True
for Nfac, Ms in [({5: 2, 13: 1, 17: 1}, [3, 4]),
                 ({13: 1, 17: 1, 29: 1, 37: 1}, [3, 4]),
                 ({5: 1, 13: 1, 17: 1, 29: 1, 37: 1}, [3, 4]),
                 ({5: 3, 401: 1, 409: 1}, [3, 4]),
                 ({313: 1, 317: 1, 337: 1, 349: 1}, [3, 4, 5]),
                 ({13: 2, 17: 2, 101: 1}, [3, 4])]:
    st = run(Nfac, Ms)
    good = st["minratio_char"] >= -1e-9 and st["minratio_mon"] >= -1e-9 and st["k_int_err"] < 1e-8
    allok &= good
    print(Nfac, Ms, {k: (float(f"{v:.3g}") if isinstance(v, float) else v) for k, v in st.items()}, "OK" if good else "VIOLATION")
print("ALL ACTUAL-CLUSTER MEMBERSHIP CHECKS PASSED" if allok else "FAILURE")
