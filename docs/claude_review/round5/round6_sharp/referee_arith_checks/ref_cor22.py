"""Referee check of sharp.md Cor. 2.2 at the level of actual residues (no dlog shortcuts in the
brute-force part).  Independent code (does not import sharp_checks).

Part A (exhaustive).  Small profiles (M = 2..4 rows, exponents e_j in {1,2}, windings k_x in Z/4,
actual primes p_j = 1 mod 4 with mixed Legendre symbols), small Q = q^a.  Enumerate EVERY tuple
(r_j) with Norm r_j = p_j mod Q (one representative of each +-pair; -r_j changes every rho_x by
the same sign (-1)^(e_j)), form rho_x = i^(-k_x) prod_j r_j^(a_xj) conj(r_j)^(e_j - a_xj) mod Q,
and record the class-difference vector (log_g(rho_x/rho_0))_(x>=1) in (Z/m)^(M-1).
  * twins profiles: claim  realisable set == {d : d_x = sigma(x) - sigma(0) mod 2};
  * profiles without twins: claim  realisable set  is contained in that parity set (necessity);
    we also record whether it is all of it.
Part B (constructive, explicit preimage (2.2)).  Random twin profiles (M = 6..14, e_j in {1,2},
random windings k), all odd prime powers Q <= 2500 with q not dividing N: random parity-respecting
target kappa, ell from (2.2), r_j by Hilbert 90 + Hensel square root (own code), then check
Norm r_j = p_j, rho_x == B g^(kappa_x - kappa_0) exactly with B = prod conj(r_j)^(e_j), and the same
at every lower level with the reduced generator (level compatibility)."""
import itertools, random, sys
import numpy as np

rng = random.Random(12345)


def is_prime(n):
    if n < 2:
        return False
    i = 2
    while i * i <= n:
        if n % i == 0:
            return False
        i += 1
    return True


def gm(u, v, Q):
    return ((u[0] * v[0] - u[1] * v[1]) % Q, (u[0] * v[1] + u[1] * v[0]) % Q)


def gc(u, Q):
    return (u[0] % Q, (-u[1]) % Q)


def gp(u, e, Q):
    r = (1, 0)
    b = (u[0] % Q, u[1] % Q)
    while e:
        if e & 1:
            r = gm(r, b, Q)
        b = gm(b, b, Q)
        e >>= 1
    return r


def gnorm(u, Q):
    return (u[0] * u[0] + u[1] * u[1]) % Q


def ginv(u, Q):
    n = pow(gnorm(u, Q), -1, Q)
    c = gc(u, Q)
    return ((c[0] * n) % Q, (c[1] * n) % Q)


def leg(n, q):
    n %= q
    if n == 0:
        return 0
    return 1 if pow(n, (q - 1) // 2, q) == 1 else -1


def mT(q, a):
    return (q - (1 if q % 4 == 1 else -1)) * q ** (a - 1)


def primefactors(n):
    out, d = [], 2
    while d * d <= n:
        if n % d == 0:
            out.append(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        out.append(n)
    return out


def generator(q, a):
    Q, m = q ** a, mT(q, a)
    fs = primefactors(m)
    for x in range(Q):
        for y in range(Q):
            if (x * x + y * y) % Q == 1:
                if all(gp((x, y), m // p, Q) != (1, 0) for p in fs):
                    return (x, y)


def dlog_table(g, m, Q):
    tab, t = {}, (1, 0)
    for k in range(m):
        tab[t] = k
        t = gm(t, g, Q)
    assert t == (1, 0)
    return tab


ps1 = [p for p in range(5, 4000) if p % 4 == 1 and is_prime(p)]

# ---------------------------------------------------------------- Part A
def twin_profile(M, extra, e2):
    """rows x 0..M-1; columns: for each x a base column b (random 0..e) and twins f(x), s(x)
    differing only in row x by eps_x; plus `extra` random columns.  e2: allow exponent 2."""
    cols, ecol, f, s, eps = [], [], {}, {}, {}
    for x in range(M):
        e = rng.choice([1, 2]) if e2 else 1
        base = [rng.randint(0, e) for _ in range(M)]
        ep = rng.choice([1, -1])
        # f column = base with row x moved by ep (stay within [0, e])
        if base[x] + ep < 0 or base[x] + ep > e:
            ep = -ep
        colf = list(base)
        colf[x] += ep
        f[x], s[x], eps[x] = len(cols), len(cols) + 1, ep
        cols += [colf, list(base)]
        ecol += [e, e]
    for _ in range(extra):
        e = rng.choice([1, 2]) if e2 else 1
        cols.append([rng.randint(0, e) for _ in range(M)])
        ecol.append(e)
    Amat = np.array(cols, dtype=np.int64).T
    return Amat, ecol, f, s, eps


def notwin_profile(M, r, e2):
    cols, ecol = [], []
    for _ in range(r):
        e = rng.choice([1, 2]) if e2 else 1
        cols.append([rng.randint(0, e) for _ in range(M)])
        ecol.append(e)
    return np.array(cols, dtype=np.int64).T, ecol


def partA_case(q, a, Amat, ecol, kvec, plist):
    Q, m = q ** a, mT(q, a)
    g = generator(q, a)
    tab = dlog_table(g, m, Q)
    M, r = Amat.shape
    # representatives (mod +-) of {r : Norm r = p_j}
    reps = []
    for j in range(r):
        allr = [(x, y) for x in range(Q) for y in range(Q) if (x * x + y * y - plist[j]) % Q == 0]
        assert len(allr) == m, (Q, len(allr), m)
        seen, rep = set(), []
        for u in allr:
            if u in seen:
                continue
            seen.add(u)
            seen.add(((-u[0]) % Q, (-u[1]) % Q))
            rep.append(u)
        reps.append(rep)
    # precompute per column j and row x the factor r^(a) conj(r)^(e-a)
    fac = [[[gm(gp(u, int(Amat[x, j]), Q), gp(gc(u, Q), ecol[j] - int(Amat[x, j]), Q), Q)
             for u in reps[j]] for x in range(M)] for j in range(r)]
    ipow = [gp((0, 1), (-kvec[x]) % 4, Q) for x in range(M)]
    found = set()
    for choice in itertools.product(*[range(len(reps[j])) for j in range(r)]):
        rho = []
        for x in range(M):
            u = ipow[x]
            for j in range(r):
                u = gm(u, fac[j][x][choice[j]], Q)
            rho.append(u)
        i0 = ginv(rho[0], Q)
        found.add(tuple(tab[gm(rho[x], i0, Q)] for x in range(1, M)))
    nu = [1 if leg(p, q) == -1 else 0 for p in plist]
    lam = 1 if q % 8 in (3, 5) else 0
    sig = [(int(sum(Amat[x, j] * nu[j] for j in range(r))) + lam * kvec[x]) % 2 for x in range(M)]
    par = set(itertools.product(*[[d for d in range(m) if (d - sig[x] + sig[0]) % 2 == 0]
                                  for x in range(1, M)]))
    return found, par


okA = True
logA = []
casesA = [(3, 1), (5, 1), (7, 1), (9, 2), (11, 1), (13, 1), (25, 2), (27, 3), (49, 2)]
for (Qv, a) in casesA:
    q = round(Qv ** (1 / a))
    m = mT(q, a)
    for trial in range(4):
        M = 3 if (m // 2) ** 6 <= 3 * 10 ** 5 else 2
        if m <= 4 and trial >= 2:
            M = 4
        Amat, ecol, f, s, eps = twin_profile(M, 0, e2=(trial % 2 == 1))
        r = Amat.shape[1]
        if (m // 2) ** r > 3 * 10 ** 6:
            continue
        plist = rng.sample([p for p in ps1 if p % q], r)
        kvec = [rng.randrange(4) for _ in range(M)]
        found, par = partA_case(q, a, Amat, ecol, kvec, plist)
        ok = (found == par)
        okA &= ok
        logA.append(("twins", Qv, M, r, ecol, kvec, len(found), len(par), ok))
        # no twins
        Amat2, ecol2 = notwin_profile(M, min(r, 4), e2=(trial % 2 == 1))
        r2 = Amat2.shape[1]
        if (m // 2) ** r2 <= 3 * 10 ** 6:
            plist2 = rng.sample([p for p in ps1 if p % q], r2)
            found2, par2 = partA_case(q, a, Amat2, ecol2, kvec, plist2)
            okA &= found2 <= par2
            logA.append(("no twins", Qv, M, r2, ecol2, kvec, len(found2), len(par2), found2 <= par2))

print("Part A (exhaustive over all residue tuples with the prescribed norms):")
for row in logA:
    print("  ", row)
print("  Part A:", "PASSED" if okA else "FAILED")

# ---------------------------------------------------------------- Part B
def hilbert90(t, q, Q):
    while True:
        u = (rng.randrange(Q), rng.randrange(Q))
        rr = ((u[0] + gm(t, gc(u, Q), Q)[0]) % Q, (u[1] + gm(t, gc(u, Q), Q)[1]) % Q)
        if gnorm(rr, q) % q:
            return rr


def sqrt_mod(n, q, a):
    n %= q ** a
    s = next(x for x in range(1, q) if (x * x - n) % q == 0)
    mod = q
    for _ in range(1, a):
        mod *= q
        s = (s - (s * s - n) * pow(2 * s, -1, mod)) % mod
    assert (s * s - n) % q ** a == 0
    return s


okB = True
nB = 0
nlev = 0
oddprimes = [q for q in range(3, 2501) if is_prime(q)]
for trial in range(6):
    M = rng.choice([6, 8, 10, 12, 14])
    Amat, ecol, f, s, eps = twin_profile(M, extra=M, e2=True)
    r = Amat.shape[1]
    plist = rng.sample(ps1[50:], r)
    kvec = [rng.randrange(4) for _ in range(M)]
    for q in oddprimes:
        if any(p % q == 0 for p in plist):
            continue
        a = 1
        while q ** (a + 1) <= 2500:
            a += 1
        if q > 400 and trial >= 2:
            continue   # keep the run short; all q <= 2500 for the first two profiles
        Q, m = q ** a, mT(q, a)
        g = generator(q, a) if Q < 400 else None
        if g is None:
            # random search for a generator (big Q)
            fs = primefactors(m)
            while True:
                u = (rng.randrange(Q), rng.randrange(Q))
                if gnorm(u, q) % q == 0:
                    continue
                t = gm(u, ginv(gc(u, Q), Q), Q)
                if all(gp(t, m // p, Q) != (1, 0) for p in fs):
                    g = t
                    break
        # iota = log_g i  (search by baby-step on the 4 candidates m/4 * {1,3})
        iota = next(c for c in (m // 4, 3 * m // 4) if gp(g, c, Q) == (0, 1))
        nu = [1 if leg(p, q) == -1 else 0 for p in plist]
        lam = 1 if q % 8 in (3, 5) else 0
        assert iota % 2 == lam
        sig = [(int(sum(Amat[x, j] * nu[j] for j in range(r))) + lam * kvec[x]) % 2 for x in range(M)]
        c0 = rng.randrange(2)
        kappa = [(sig[x] + c0 + 2 * rng.randrange(m // 2)) % m for x in range(M)]
        kappa0 = (c0 + 2 * rng.randrange(m // 2)) % m          # any kappa_0 of the right parity
        Anu = [int(sum(Amat[x, j] * nu[j] for j in range(r))) for x in range(M)]
        y = []
        for x in range(M):
            num = (kappa[x] - kappa0 + iota * kvec[x] - Anu[x]) % m
            okB &= num % 2 == 0
            y.append(num // 2)
        ellp = [0] * r
        for x in range(M):
            ellp[f[x]] += y[x] * eps[x]
            ellp[s[x]] -= y[x] * eps[x]
        ell = [(nu[j] + 2 * ellp[j]) % m for j in range(r)]
        okB &= all(ell[j] % 2 == nu[j] for j in range(r))
        rj = []
        for j in range(r):
            t = gp(g, ell[j], Q)
            r0 = hilbert90(t, q, Q)
            sc = sqrt_mod(plist[j] * pow(gnorm(r0, Q), -1, Q), q, a)
            rr = ((r0[0] * sc) % Q, (r0[1] * sc) % Q)
            okB &= gnorm(rr, Q) == plist[j] % Q and gm(rr, ginv(gc(rr, Q), Q), Q) == t
            rj.append(rr)
        B = (1, 0)
        for j in range(r):
            B = gm(B, gp(gc(rj[j], Q), ecol[j], Q), Q)
        for lev in range(1, a + 1):
            Ql, ml = q ** lev, mT(q, lev)
            gl = (g[0] % Ql, g[1] % Ql)
            Bl = (B[0] % Ql, B[1] % Ql)
            for x in range(M):
                u = gp((0, 1), (-kvec[x]) % 4, Ql)
                for j in range(r):
                    u = gm(u, gm(gp(rj[j], int(Amat[x, j]), Ql), gp(gc(rj[j], Ql), ecol[j] - int(Amat[x, j]), Ql), Ql), Ql)
                pred = gm(Bl, gp(gl, (kappa[x] - kappa0) % ml, Ql), Ql)
                okB &= (u == pred)
            nlev += 1
        nB += 1
print(f"Part B: {nB} (profile, prime) cases, {nlev} levels checked:", "PASSED" if okB else "FAILED")
print("ALL COR 2.2 CHECKS PASSED" if (okA and okB) else "FAILURE")
