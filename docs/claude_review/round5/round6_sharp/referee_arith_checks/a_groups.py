"""Referee check (arith lens) of sharp.md Prop. 2.1, items 1-6, EXHAUSTIVELY over Z[i]/Q for every
odd prime power Q = q^a <= QMAX (default 1000).  Independent numpy code (does not import gres.py).

For each Q:
 (1) T_Q = {u : Norm u = 1 mod Q} has m = (q - chi(q)) q^(a-1) elements, 4 | m, and is cyclic
     (an element of order m exists and its powers exhaust T_Q);
 (2) r -> r/conj(r) (= r^2/Norm r) maps G_Q = (Z[i]/Q)^* ONTO T_Q, every fibre has phi(Q)
     elements, and the fibre of t is exactly r_0 (Z/Q)^* (checked for every t);
 (3) Legendre(Norm r) is constant on every fibre; it equals -1 iff log_g t is odd (t not in T^2);
 (4) for EVERY t in T_Q and EVERY n in (Z/Q)^*: #{r : r/conj r = t, Norm r = n} is 2 if
     nu(t) = (n/q) and 0 otherwise  (so: exists iff nu(t) = (n/q));
 (5) i in T^2 iff q = +-1 mod 8;  -1 in T^2 always;
 (6) (a >= 2) reduction T_{q^a} -> T_{q^(a-1)} is onto, maps a generator to a generator, and
     nu_{q^a}(t) = nu_{q^(a-1)}(red t) with INDEPENDENTLY chosen generators at the two levels.
"""
import sys, math, random
import numpy as np

QMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 1000


def primes_upto(n):
    s = np.ones(n + 1, dtype=bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return [int(i) for i in np.nonzero(s)[0]]


def factor(n):
    f, d = {}, 2
    while d * d <= n:
        while n % d == 0:
            f[d] = f.get(d, 0) + 1; n //= d
        d += 1
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


def gmul(u, v, Q):
    return ((u[0] * v[0] - u[1] * v[1]) % Q, (u[0] * v[1] + u[1] * v[0]) % Q)


def gpow(u, e, Q):
    r, b = (1, 0), (u[0] % Q, u[1] % Q)
    while e:
        if e & 1:
            r = gmul(r, b, Q)
        b = gmul(b, b, Q); e >>= 1
    return r


class Level:
    """all data of Z[i]/Q, Q = q^a"""

    def __init__(self, q, a, seed):
        self.q, self.a = q, a
        Q = q ** a; self.Q = Q
        chi = 1 if q % 4 == 1 else -1
        self.m = (q - chi) * q ** (a - 1)
        x = np.repeat(np.arange(Q, dtype=np.int64), Q)
        y = np.tile(np.arange(Q, dtype=np.int64), Q)
        self.x, self.y = x, y
        nrm = (x * x + y * y) % Q
        self.nrm = nrm
        self.unit = (nrm % q) != 0
        # inverse table on Z/Q
        inv = np.zeros(Q, dtype=np.int64)
        for n in range(1, Q):
            if n % q:
                inv[n] = pow(n, -1, Q)
        self.inv = inv
        leg = np.zeros(Q, dtype=np.int64)
        for n in range(1, Q):
            if n % q:
                leg[n] = 1 if pow(n % q, (q - 1) // 2, q) == 1 else -1
        self.leg = leg
        Tmask = nrm == 1
        self.Tidx = np.nonzero(Tmask)[0]
        # generator
        rng = random.Random(seed)
        fm = factor(self.m)
        Tlist = self.Tidx.tolist()
        while True:
            t = rng.choice(Tlist)
            g = (t // Q, t % Q)
            if all(gpow(g, self.m // p, Q) != (1, 0) for p in fm):
                break
        self.g = g
        dl = -np.ones(Q * Q, dtype=np.int64)
        t = (1, 0)
        for k in range(self.m):
            idx = t[0] * Q + t[1]
            if dl[idx] != -1:
                raise RuntimeError("generator has small order")
            dl[idx] = k
            t = gmul(t, g, Q)
        assert t == (1, 0)
        self.dl = dl
        # phi(r) = r^2 / Norm(r) for units
        u = self.unit
        xu, yu = x[u], y[u]
        ni = inv[nrm[u]]
        px = ((xu * xu - yu * yu) % Q) * ni % Q
        py = ((2 * xu * yu) % Q) * ni % Q
        self.phi_u = px * Q + py          # index of r/conj r, for units in order
        self.uidx = np.nonzero(u)[0]


def check(q, a, report):
    L = Level(q, a, seed=1000 * q + a)
    Q, m = L.Q, L.m
    phiQ = Q - Q // q
    ok = {}
    # (1)
    ok['1_order'] = len(L.Tidx) == m and m % 4 == 0
    ok['1_cyclic'] = bool(np.all(L.dl[L.Tidx] >= 0)) and int((L.dl >= 0).sum()) == m
    # (2) onto + fibre sizes
    img = L.dl[L.phi_u]
    ok['2_in_T'] = bool(np.all(img >= 0))
    cnt = np.bincount(img, minlength=m)
    ok['2_onto_fibres'] = bool(np.all(cnt == phiQ))
    # fibre = r0 (Z/Q)^*  : for every t pick r0 and check phi(u r0) = t for all rational units u,
    # then sizes force equality
    first = {}
    order = np.argsort(img, kind='stable')
    starts = np.searchsorted(img[order], np.arange(m))
    r0idx = L.uidx[order[starts]]
    rats = np.array([u for u in range(1, Q) if u % q], dtype=np.int64)
    X0, Y0 = r0idx // Q, r0idx % Q
    good = True
    for u in rats:
        xx, yy = (X0 * u) % Q, (Y0 * u) % Q
        nn = (xx * xx + yy * yy) % Q
        ni = L.inv[nn]
        px = ((xx * xx - yy * yy) % Q) * ni % Q
        py = ((2 * xx * yy) % Q) * ni % Q
        good &= bool(np.all(L.dl[px * Q + py] == np.arange(m)))
    ok['2_fibre_is_coset'] = good
    # (3) nu constant on fibres, = parity of dlog
    legs = L.leg[L.nrm[L.uidx]]
    mn = np.full(m, 2); mx = np.full(m, -2)
    np.minimum.at(mn, img, legs); np.maximum.at(mx, img, legs)
    ok['3_nu_well_defined'] = bool(np.all(mn == mx))
    nu = mn
    ok['3_nu_is_square_class'] = bool(np.all((nu == 1) == (np.arange(m) % 2 == 0)))
    ok['3_onto'] = set(nu.tolist()) == {1, -1}
    # (4) all t, all n
    H = np.zeros((m, Q), dtype=np.int64)
    np.add.at(H, (img, L.nrm[L.uidx]), 1)
    nvals = np.array([n for n in range(Q) if n % q], dtype=np.int64)
    Hs = H[:, nvals]
    pred = (nu[:, None] == L.leg[nvals][None, :])
    ok['4_exists_iff'] = bool(np.all((Hs > 0) == pred))
    ok['4_count_is_2'] = bool(np.all(Hs[pred] == 2))
    # (5)
    di = int(L.dl[0 * Q + 1]); dm1 = int(L.dl[(Q - 1) * Q + 0])
    ok['5_i'] = (di % 2 == 0) == (q % 8 in (1, 7))
    ok['5_minus1'] = dm1 % 2 == 0 and dm1 == m // 2
    # (6)
    if a >= 2:
        Lo = Level(q, a - 1, seed=777 * q + a)   # independent generator
        Qo = Lo.Q
        Tx, Ty = L.Tidx // Q, L.Tidx % Q
        red = (Tx % Qo) * Qo + (Ty % Qo)
        ok['6_into'] = bool(np.all(Lo.dl[red] >= 0))
        ok['6_onto'] = len(set(Lo.dl[red].tolist())) == Lo.m
        gred = (L.g[0] % Qo, L.g[1] % Qo)
        o = Lo.m
        for p in factor(Lo.m):
            while o % p == 0 and gpow(gred, o // p, Qo) == (1, 0):
                o //= p
        ok['6_gen_to_gen'] = o == Lo.m
        ok['6_nu_compat'] = bool(np.all((L.dl[L.Tidx] % 2) == (Lo.dl[red] % 2)))
        # with compatible generator: dlog reduces mod m_lower
        gi = Lo.dl[gred[0] * Qo + gred[1]]          # log of red(g) in lower base
        ok['6_dlog_reduces'] = bool(np.all(Lo.dl[red] == (L.dl[L.Tidx] * gi) % Lo.m))
    report.append((Q, q, a, m, q % 8, ok))
    return all(ok.values())


if __name__ == "__main__":
    allok = True
    report = []
    mods = []
    for q in primes_upto(QMAX):
        if q == 2:
            continue
        a = 1
        while q ** a <= QMAX:
            mods.append((q, a)); a += 1
    for (q, a) in mods:
        allok &= check(q, a, report)
    keys = sorted({k for r in report for k in r[5]})
    fails = [(r[0], [k for k, v in r[5].items() if not v]) for r in report if not all(r[5].values())]
    print(f"odd prime powers Q <= {QMAX}: {len(mods)}  (prime powers with a >= 2: {sum(1 for _, a in mods if a >= 2)})")
    print("items checked:", keys)
    print("failures:", fails if fails else "none")
    # count moduli <= 400 for comparison with sharp.md's '89'
    print("odd prime powers <= 400:", sum(1 for (q, a) in mods if q ** a <= 400))
    print("ALL PROP 2.1 CHECKS PASSED" if allok else "PROP 2.1 FAILURE")
