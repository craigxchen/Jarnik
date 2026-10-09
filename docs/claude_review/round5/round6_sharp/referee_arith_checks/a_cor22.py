"""Referee check (arith lens) of sharp.md Cor. 2.2 (items 1-5) and Remarks 2.3, independent code.

Profiles (general exponents allowed):  P1: e = 1, M = 6, twins + random columns;
P2: exponents e_j in {1,2,3}, M = 5, twins (with e_j up to 2) + random columns.
Primes: actual p_j = 1 mod 4 (> 10^4), pi_j from sum of two squares.
Moduli: every odd prime q <= X (X = 130), top power Q = q^(A_q), all levels a <= A_q.

Checks:
 C1  (item 1) for each j: {log_g(r/conj r) : Norm r = p_j mod Q} = {ell : ell = nu_j mod 2},
     each value attained by exactly two r (r, -r); pi_j mod Q has the forced parity;
     nu_j(q) = [(p_j/q) = -1] = [(q/p_j) = -1].
 C2  (items 2-3) random admissible r_j, random windings k (|k| <= 9): rho_x = B g^kappa_x with B
     common and kappa = A ell - iota k; for EVERY c in {-2..2}^M with sum c = 0, every level a,
     every s: prod rho^c = i^s mod q^a  <=>  sum c kappa = iota s mod m_(q^a).
 C3  (item 4) kappa_x - kappa_y = sigma_q(x) - sigma_q(y) mod 2, sigma with lambda_q k_x term,
     lambda_q = [q = +-3 mod 8] = iota mod 2.
 C4  (item 5) random parity-respecting kappa* (random constant offset), random windings: the
     explicit preimage (2.2) gives ell = nu mod 2, residues r_j with Norm r_j = p_j mod Q exist
     (Hilbert 90 + square root, built here), and the realised point classes are kappa* - kappa_0.
 C5  exhaustive at tiny size (M = 3, 7 columns, q in {3,5,7,13,17}): the set of realisable class
     vectors {A ell - iota k : ell = nu mod 2} equals {kappa : kappa = A nu - iota k mod 2}.
 C6  Remark 2.3(2): for 8 | m, ell(pi) mod 4 is invariant under change of generator, of
     associate and of conjugate (when (p/q) = 1); it is not a function of (p/q); and for inert
     q = 7 mod 8 it equals the quartic symbol of -q modulo pi (quartic reciprocity), i.e. it is
     determined by data at the modulus p, not at q.
"""
import random, math, itertools, sys

rng = random.Random(20261009)


def is_prime(n):
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


def primes_upto(n):
    return [p for p in range(2, n + 1) if is_prime(p)]


def two_squares(p):
    for x in range(1, int(math.isqrt(p)) + 1):
        y2 = p - x * x
        y = math.isqrt(y2)
        if y * y == y2:
            return (x, y)


def leg(a, q):
    a %= q
    if a == 0:
        return 0
    return 1 if pow(a, (q - 1) // 2, q) == 1 else -1


def gm(u, v, Q):
    return ((u[0] * v[0] - u[1] * v[1]) % Q, (u[0] * v[1] + u[1] * v[0]) % Q)


def gp(u, e, Q):
    if e < 0:
        u = ginv(u, Q); e = -e
    r, b = (1, 0), (u[0] % Q, u[1] % Q)
    while e:
        if e & 1:
            r = gm(r, b, Q)
        b = gm(b, b, Q); e >>= 1
    return r


def gconj(u, Q):
    return (u[0] % Q, (-u[1]) % Q)


def gnorm(u, Q):
    return (u[0] * u[0] + u[1] * u[1]) % Q


def ginv(u, Q):
    n = pow(gnorm(u, Q), -1, Q)
    c = gconj(u, Q)
    return (c[0] * n % Q, c[1] * n % Q)


def factor(n):
    f, d = {}, 2
    while d * d <= n:
        while n % d == 0:
            f[d] = f.get(d, 0) + 1; n //= d
        d += 1
    if n > 1:
        f[n] = 1
    return f


class Mod:
    def __init__(self, q, A):
        self.q, self.A = q, A
        self.Q = q ** A
        chi = 1 if q % 4 == 1 else -1
        self.m = (q - chi) * q ** (A - 1)
        Q, m = self.Q, self.m
        T = [(x, y) for x in range(Q) for y in range(Q) if (x * x + y * y) % Q == 1]
        assert len(T) == m
        fm = factor(m)
        while True:
            g = rng.choice(T)
            if all(gp(g, m // p, Q) != (1, 0) for p in fm):
                break
        self.g = g
        self.dl = {}
        t = (1, 0)
        for k in range(m):
            self.dl[t] = k
            t = gm(t, g, Q)
        self.iota = self.dl[(0, 1)]

    def m_at(self, a):
        chi = 1 if self.q % 4 == 1 else -1
        return (self.q - chi) * self.q ** (a - 1)

    def t_of(self, r):
        return gm(r, ginv(gconj(r, self.Q), self.Q), self.Q)


def sqrt_mod(n, q, A):
    Q = q ** A
    n %= Q
    s = next((x for x in range(1, q) if (x * x - n) % q == 0), None)
    if s is None:
        return None
    mod = q
    for _ in range(1, A):
        mod *= q
        s = (s - (s * s - n) * pow(2 * s, -1, mod)) % mod
    assert (s * s - n) % Q == 0
    return s


def residue_with(t, n, md):
    """r with r/conj r = t and Norm r = n mod Q (Hilbert 90 then rescale)"""
    Q, q = md.Q, md.q
    while True:
        u = (rng.randrange(Q), rng.randrange(Q))
        tu = gm(t, gconj(u, Q), Q)
        r0 = ((u[0] + tu[0]) % Q, (u[1] + tu[1]) % Q)
        if gnorm(r0, Q) % q:
            break
    assert md.t_of(r0) == t
    s = sqrt_mod(n * pow(gnorm(r0, Q), -1, Q), q, md.A)
    if s is None:
        return None
    r = (r0[0] * s % Q, r0[1] * s % Q)
    assert gnorm(r, Q) == n % Q and md.t_of(r) == t
    return r


def build_profile(M, extra, emax_extra, emax_twin):
    """rows a_x; returns (rows, e, twins[(f, s, eps)])"""
    while True:
        cols, e, twins = [], [], []
        for x in range(M):
            ej = rng.randint(1, emax_twin)
            s = [rng.randint(0, ej) for _ in range(M)]
            if s[x] < ej:
                eps = 1
            else:
                eps = -1
            f = list(s); f[x] += eps
            twins.append((len(cols), len(cols) + 1, eps))
            cols.append(f); cols.append(s); e += [ej, ej]
        for _ in range(extra):
            ej = rng.randint(1, emax_extra)
            cols.append([rng.randint(0, ej) for _ in range(M)]); e.append(ej)
        # normalisation: min 0, max e_j for each column; rows distinct
        okn = all(min(c) == 0 and max(c) == ej for c, ej in zip(cols, e))
        rows = [tuple(c[x] for c in cols) for x in range(M)]
        if okn and len(set(rows)) == M:
            return [list(r) for r in rows], e, twins


def run_profile(name, M, extra, emax_extra, emax_twin, X, ncheck_c=None):
    rows, e, twins = build_profile(M, extra, emax_extra, emax_twin)
    r = len(e)
    ps, n = [], 10 ** 4 + 1
    while len(ps) < r:
        if n % 4 == 1 and is_prime(n):
            ps.append(n)
        n += 4 * rng.randint(1, 30)
    pis = [two_squares(p) for p in ps]
    ok = {k: True for k in ["C1", "C2", "C3", "C4"]}
    nmod = ncols = 0
    cs = [c for c in itertools.product(range(-2, 3), repeat=M) if sum(c) == 0 and any(c)]
    for q in primes_upto(X):
        if q == 2:
            continue
        A = 1
        while q ** (A + 1) <= X:
            A += 1
        md = Mod(q, A)
        Q, m, iota = md.Q, md.m, md.iota
        nmod += 1
        nu = [1 if leg(p, q) == -1 else 0 for p in ps]
        ok["C1"] &= all(nu[j] == (1 if leg(q, ps[j]) == -1 else 0) for j in range(r))
        lam = 1 if q % 8 in (3, 5) else 0
        ok["C3"] &= (iota % 2 == lam)
        # C1: attainable ell for each j
        normset = {}
        for x in range(Q):
            for y in range(Q):
                nn = (x * x + y * y) % Q
                if nn % q:
                    normset.setdefault(nn, []).append((x, y))
        choices = []
        for j, p in enumerate(ps):
            rs = normset[p % Q]
            ells = [md.dl[md.t_of(rr)] for rr in rs]
            cnt = {}
            for l in ells:
                cnt[l] = cnt.get(l, 0) + 1
            ok["C1"] &= set(cnt) == {l for l in range(m) if l % 2 == nu[j]}
            ok["C1"] &= all(v == 2 for v in cnt.values())
            pij = (pis[j][0] % Q, pis[j][1] % Q)
            ok["C1"] &= md.dl[md.t_of(pij)] % 2 == nu[j]
            choices.append(rs)
        # C2/C3: random admissible residues and windings
        for trial in range(2):
            rj = [rng.choice(choices[j]) for j in range(r)]
            ell = [md.dl[md.t_of(x)] for x in rj]
            k = [rng.randint(-9, 9) for _ in range(M)]
            rho = []
            for x in range(M):
                v = gp((0, 1), -k[x], Q)
                for j in range(r):
                    v = gm(v, gm(gp(rj[j], rows[x][j], Q), gp(gconj(rj[j], Q), e[j] - rows[x][j], Q), Q), Q)
                rho.append(v)
            kap = [(sum(rows[x][j] * ell[j] for j in range(r)) - iota * k[x]) % m for x in range(M)]
            Bs = {gm(rho[x], gp(md.g, -kap[x], Q), Q) for x in range(M)}
            ok["C2"] &= len(Bs) == 1
            sig = [(sum(rows[x][j] * nu[j] for j in range(r)) + lam * k[x]) % 2 for x in range(M)]
            ok["C3"] &= all((kap[x] - kap[0] - sig[x] + sig[0]) % 2 == 0 for x in range(M))
            # collisions at every level, every unit
            sub = cs if ncheck_c is None else rng.sample(cs, min(ncheck_c, len(cs)))
            for a in range(1, A + 1):
                qa, ma = q ** a, md.m_at(a)
                rho_a = [(v[0] % qa, v[1] % qa) for v in rho]
                units = [gp((0, 1), s, qa) for s in range(4)]
                for c in sub:
                    P = (1, 0)
                    for x in range(M):
                        if c[x]:
                            P = gm(P, gp(rho_a[x], c[x], qa), qa)
                    sk = sum(c[x] * kap[x] for x in range(M))
                    for s in range(4):
                        ok["C2"] &= (P == units[s]) == ((sk - iota * s) % ma == 0)
                    ncols += 1
        # C4: explicit preimage (2.2)
        for trial in range(3):
            k = [rng.randint(-9, 9) for _ in range(M)]
            sig = [(sum(rows[x][j] * nu[j] for j in range(r)) + lam * k[x]) % 2 for x in range(M)]
            off = rng.randint(0, 1)
            kst = [(sig[x] + off + 2 * rng.randrange(m // 2)) % m for x in range(M)]
            kap0 = off  # kappa* - kappa_0 + iota k - A nu must be even
            Anu = [sum(rows[x][j] * nu[j] for j in range(r)) for x in range(M)]
            ydiff = [(kst[x] - kap0 + iota * k[x] - Anu[x]) for x in range(M)]
            ok["C4"] &= all(v % 2 == 0 for v in ydiff)
            y = [(v // 2) % (m // 2) for v in ydiff]
            ellp = [0] * r
            for x, (f, s_, eps) in enumerate(twins):
                ellp[f] += y[x] * eps
                ellp[s_] -= y[x] * eps
            ell = [(nu[j] + 2 * ellp[j]) % m for j in range(r)]
            ok["C4"] &= all(ell[j] % 2 == nu[j] for j in range(r))
            rj = []
            for j in range(r):
                rr = residue_with(gp(md.g, ell[j], Q), ps[j], md)
                ok["C4"] &= rr is not None
                rj.append(rr)
            rho = []
            for x in range(M):
                v = gp((0, 1), -k[x], Q)
                for j in range(r):
                    v = gm(v, gm(gp(rj[j], rows[x][j], Q), gp(gconj(rj[j], Q), e[j] - rows[x][j], Q), Q), Q)
                rho.append(v)
            B = gm(rho[0], gp(md.g, -(kst[0] - kap0), Q), Q)
            ok["C4"] &= all(gm(rho[x], gp(md.g, -(kst[x] - kap0), Q), Q) == B for x in range(M))
            ok["C4"] &= all(gnorm(rj[j], Q) == ps[j] % Q for j in range(r))
    print(f"{name}: M={M} r={r} e={e}  moduli={nmod}  character-level tests={ncols}  ->",
          {k: v for k, v in ok.items()})
    return all(ok.values())


def tiny_exhaustive():
    good = True
    for q in (3, 5, 7, 13, 17):
        md = Mod(q, 1)
        m, iota = md.m, md.iota
        M = 3
        rows, e, twins = build_profile(M, 1, 1, 1)
        r = len(e)
        nu = [rng.randint(0, 1) for _ in range(r)]
        k = [rng.randint(-3, 3) for _ in range(M)]
        S = set()
        for lp in itertools.product(range(m // 2), repeat=r):
            ell = [nu[j] + 2 * lp[j] for j in range(r)]
            S.add(tuple((sum(rows[x][j] * ell[j] for j in range(r)) - iota * k[x]) % m for x in range(M)))
        target = {kv for kv in itertools.product(range(m), repeat=M)
                  if all((kv[x] - sum(rows[x][j] * nu[j] for j in range(r)) + iota * k[x]) % 2 == 0 for x in range(M))}
        good &= (S == target)
        print(f"  C5 q={q} m={m} r={r}: |realisable|={len(S)}  |parity-respecting|={len(target)}  equal={S == target}")
    return good


def quartic():
    good = True
    out = []
    gps = [(p, two_squares(p)) for p in primes_upto(3000) if p % 4 == 1]
    for q in [7, 17, 23, 31, 41, 47, 71, 73, 79, 89, 97, 103]:
        md = Mod(q, 1)
        Q, m = md.Q, md.m
        if m % 8:
            continue
        # second generator g2 = g^u, u odd unit mod m
        u = next(u for u in range(3, m, 2) if math.gcd(u, m) == 1 and u % 4 == 3)
        g2 = gp(md.g, u, Q)
        dl2 = {}
        t = (1, 0)
        for kk in range(m):
            dl2[t] = kk; t = gm(t, g2, Q)
        vals = {1: set(), -1: set()}
        inv_ok = True
        recip_ok = True
        for p, (x, y) in gps:
            if p == q:
                continue
            L = leg(p, q)
            pis = [(x, y), (-y, x), (-x, -y), (y, -x), (x, -y)]   # associates and a conjugate
            ells = [md.dl[md.t_of((a % Q, b % Q))] for a, b in pis]
            if L == 1:
                inv_ok &= len({l % 4 for l in ells}) == 1 and len({dl2[md.t_of((a % Q, b % Q))] % 4 for a, b in pis}) == 1
                vals[1].add(ells[0] % 4)
                if q % 4 == 3:
                    # quartic symbol of -q mod pi: (-q)^((p-1)/4) mod p, with i -> the root of x^2=-1 mod p
                    # given by pi = x + i y = 0 mod pi  =>  i = -x/y mod p
                    ival = (-x * pow(y, -1, p)) % p
                    w = pow((-q) % p, (p - 1) // 4, p)
                    # w in {1, -1, i, -i}; with (q/p) = 1 it is +-1
                    sym = 0 if w == 1 else (2 if w == p - 1 else None)
                    recip_ok &= sym is not None and (ells[0] % 4) == sym
            else:
                vals[-1].add(ells[0] % 4)
        good &= inv_ok and recip_ok and vals[1] == {0, 2}
        out.append((q, m, inv_ok, sorted(vals[1]), sorted(vals[-1]), recip_ok if q % 4 == 3 else 'n/a'))
    for row in out:
        print(f"  C6 q={row[0]:3d} m={row[1]:3d}: invariant(mod 4)={row[2]}  ell mod 4 for (p/q)=1: {row[3]}  for -1: {row[4]}  quartic reciprocity [pi/q]_4=[-q/pi]_4: {row[5]}")
    return good


if __name__ == "__main__":
    allok = True
    allok &= run_profile("P1", 6, 8, 1, 1, 130, ncheck_c=400)
    allok &= run_profile("P2", 5, 6, 3, 2, 130, ncheck_c=None)
    print("C5 (exhaustive realisable set):")
    allok &= tiny_exhaustive()
    print("C6 (quartic part):")
    allok &= quartic()
    print("ALL COR 2.2 CHECKS PASSED" if allok else "COR 2.2 FAILURE")
