"""Referee check (independent of sharp_checks/gres.py) of sharp.md Prop. 2.1, items 1-6.

For every odd prime power Q = q^a <= QMAX, exhaustively over Z[i]/Q (numpy, exact int64):
 (1) T_Q = {u : Norm u = 1} has order m = (q - chi(q)) q^(a-1), 4 | m, and is cyclic
     (an explicit g whose powers enumerate T_Q);
 (2) t(r) = r/conj(r) = r^2 Norm(r)^(-1) maps G_Q onto T_Q, and its fibres are exactly the
     cosets r (Z/Q)^*  (via a canonical coset representative);
 (3) Legendre(Norm r) is constant on fibres; the induced nu is +1 exactly on T_Q^2, index 2;
 (4) FULL realisability table: the multiset {(t(r), Norm r) : r in G_Q} equals
     {(t, n) : nu(t) = (n/q)} with every pair hit exactly twice (all t, all n, not a sample);
 (5) i in T_Q^2 iff q = +-1 mod 8;  -1 in T_Q^2 always;  also order of i is 4;
 (6) reduction mod q^(a-1) maps T_Q onto T_(Q/q), sends EVERY generator to a generator, and
     nu_Q = nu_(Q/q) o red.
Usage: python3 ref_groups.py [QMAX]   (default 700)."""
import sys
import numpy as np

QMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 700


def primes_upto(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return np.nonzero(s)[0].tolist()


def leg(n, q):
    n = np.asarray(n) % q
    out = np.ones_like(n)
    # Euler criterion, vectorised by repeated squaring
    e = (q - 1) // 2
    r = np.ones_like(n)
    b = n.copy()
    while e:
        if e & 1:
            r = (r * b) % q
        b = (b * b) % q
        e >>= 1
    out = np.where(n == 0, 0, np.where(r == 1, 1, -1))
    return out


def factor_set(n):
    out, d = set(), 2
    while d * d <= n:
        while n % d == 0:
            out.add(d)
            n //= d
        d += 1
    if n > 1:
        out.add(n)
    return out


def mul(x1, y1, x2, y2, Q):
    return (x1 * x2 - y1 * y2) % Q, (x1 * y2 + y1 * x2) % Q


def powe(x, y, e, Q):
    rx, ry = 1, 0
    bx, by = x % Q, y % Q
    while e:
        if e & 1:
            rx, ry = mul(rx, ry, bx, by, Q)
        bx, by = mul(bx, by, bx, by, Q)
        e >>= 1
    return rx, ry


ok_all = True
rows = []
mods = []
for q in primes_upto(QMAX):
    if q == 2:
        continue
    a = 1
    while q ** a <= QMAX:
        mods.append((q, a))
        a += 1

gens_cache = {}
nu_cache = {}
for (q, a) in mods:
    Q = q ** a
    chi = 1 if q % 4 == 1 else -1
    m = (q - chi) * q ** (a - 1)
    phiQ = Q - Q // q
    X, Y = np.meshgrid(np.arange(Q, dtype=np.int64), np.arange(Q, dtype=np.int64), indexing="ij")
    X = X.ravel(); Y = Y.ravel()
    Nm = (X * X + Y * Y) % Q
    unit = (Nm % q) != 0
    ok = True
    # (1)
    Tmask = Nm == 1
    Tidx = np.nonzero(Tmask)[0]
    ok &= len(Tidx) == m and m % 4 == 0
    # find a generator: element whose powers have exact order m
    fs = factor_set(m)
    g = None
    for idx in Tidx:
        x, y = int(X[idx]), int(Y[idx])
        if all(powe(x, y, m // p, Q) != (1, 0) for p in fs):
            g = (x, y)
            break
    ok &= g is not None
    # enumerate powers -> dlog table on T
    dlog = -np.ones(Q * Q, dtype=np.int64)
    cx, cy = 1, 0
    for k in range(m):
        dlog[cx * Q + cy] = k
        cx, cy = mul(cx, cy, g[0], g[1], Q)
    ok &= (cx, cy) == (1, 0)
    ok &= np.array_equal(np.sort(np.nonzero(dlog >= 0)[0]), Tidx)   # powers of g = T
    # (2) t(r) = r^2 / Norm(r)
    U = np.nonzero(unit)[0]
    ux, uy = X[U], Y[U]
    un = Nm[U]
    inv_tab = np.zeros(Q, dtype=np.int64)
    for v in range(1, Q):
        if v % q:
            inv_tab[v] = pow(v, -1, Q)
    sx, sy = mul(ux, uy, ux, uy, Q)
    tx = (sx * inv_tab[un]) % Q
    ty = (sy * inv_tab[un]) % Q
    tid = tx * Q + ty
    ok &= bool(np.all(Tmask[tid]))                      # lands in T
    ok &= len(np.unique(tid)) == m                       # onto T
    cnt = np.bincount(tid, minlength=Q * Q)[Tidx]
    ok &= bool(np.all(cnt == phiQ))                      # fibre sizes |(Z/Q)^*|
    # canonical coset representative of r (Z/Q)^*: scale so a coordinate that is a unit becomes 1
    xunit = (ux % q) != 0
    cxr = np.where(xunit, 1, (ux * inv_tab[np.where(xunit, 1, uy)]) % Q)
    cyr = np.where(xunit, (uy * inv_tab[np.where(xunit, ux, 1)]) % Q, 1)
    rep = cxr * Q + cyr
    # fibre of t == exactly one coset: map rep -> t is well defined and injective
    pairs = np.unique(np.stack([rep, tid], 1), axis=0)
    ok &= len(np.unique(pairs[:, 0])) == len(pairs) and len(np.unique(pairs[:, 1])) == len(pairs)
    # (3) nu well defined, = squareness
    L = leg(un, q)
    nu_t = np.zeros(Q * Q, dtype=np.int64)
    nu_t[tid] = L
    # well-definedness: every (t, L) pair consistent
    tl = np.unique(np.stack([tid, L], 1), axis=0)
    ok &= len(np.unique(tl[:, 0])) == len(tl)
    dl = dlog[Tidx]
    ok &= bool(np.all((nu_t[Tidx] == 1) == (dl % 2 == 0)))
    ok &= int(np.sum(nu_t[Tidx] == 1)) == m // 2
    # (4) full realisability table
    key = tid * Q + un            # (t, n)
    hit, kcounts = np.unique(key, return_counts=True)
    ok &= bool(np.all(kcounts == 2))
    ht, hn = hit // Q, hit % Q
    ok &= bool(np.all(nu_t[ht] == leg(hn, q)))
    ok &= len(hit) == m * phiQ // 2   # every compatible (t, n) is hit
    # (5) units
    di = dlog[0 * Q + 1]          # i = (0, 1)
    dm1 = dlog[(Q - 1) * Q + 0]   # -1
    ok &= (di % 2 == 0) == (q % 8 in (1, 7))
    ok &= dm1 % 2 == 0 and dm1 == m // 2 and (4 * di) % m == 0 and (2 * di) % m == m // 2 % m
    # (6) reduction
    if a >= 2:
        Qp = Q // q
        gp, dlp, nup, mp = gens_cache[(q, a - 1)]
        redid = (X[Tidx] % Qp) * Qp + (Y[Tidx] % Qp)
        ok &= len(np.unique(redid)) == mp          # onto T_(Q/q)
        # dlog compatibility for one generator, then every generator g^k (k unit mod m)
        dred = dlp[redid]
        # red(g) has dlog d0 at level a-1; red(g^k) = red(g)^k
        d0 = dlp[(g[0] % Qp) * Qp + (g[1] % Qp)]
        ok &= np.gcd(int(d0), mp) == 1             # generator -> generator
        ok &= bool(np.all(dred == (dl * d0) % mp))
        # every generator of T_Q is g^k with gcd(k, m) = 1 and reduces to red(g)^k, a generator
        ks = [k for k in range(m) if np.gcd(k, m) == 1]
        ok &= all(np.gcd((k * int(d0)) % mp, mp) == 1 for k in ks)
        ok &= bool(np.all(nup[redid] == nu_t[Tidx]))
    gens_cache[(q, a)] = (g, dlog, nu_t, m)
    ok_all &= ok
    rows.append((Q, q, a, m, q % 8, int(di % 2 == 0), ok))

print(f"moduli: {len(mods)} odd prime powers <= {QMAX}")
bad = [r for r in rows if not r[-1]]
print("failures:", bad)
print("sample rows (Q, q, a, m, q mod 8, i square, ok):")
for r in rows[:12] + rows[-5:]:
    print("  ", r)
print("ALL GROUP CHECKS PASSED" if ok_all else "FAILURE")
