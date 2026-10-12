"""Exact check of Proposition B.3 (parity-refined pigeonhole) on ACTUAL lattice points.

For squarefree N = p_1 ... p_r (p_j = 1 mod 4) every lattice point of norm N is
    z = i^e prod_j pi_j^(a_j) conj(pi_j)^(1 - a_j),   a in {0,1}^r, e in Z/4.
For an odd prime q not dividing N put
    sigma_q(z) = sum_j a_j [(p_j/q) = -1] + e [(2/q) = -1]  mod 2.
Claims checked (integer arithmetic):
  (1) t = z_x conj(z_y) / N  mod q lies in the norm-one group T_q, and t^(m_q/2) = (-1)^(sigma(x) - sigma(y))
      (Euler criterion in the cyclic group T_q of order m_q = q - (-1/q));
  (2) hence z_x = z_y mod q only if sigma(x) = sigma(y);
  (3) #{x<y : z_x = z_y mod q^a} >= E(n_0, m/2) + E(n_1, m/2)  for every q^a <= Xmax, where
      m = m_(q^a) = (q - (-1/q)) q^(a-1), n_i = #{x : sigma(x) = i};
  (4) E(n_0, Q) + E(n_1, Q) >= E(n_0 + n_1, 2Q) and <= E(n_0 + n_1, Q) (exhaustive small table).
Clusters: random subsets and arcs (consecutive arguments) of the full point set.
"""
import math, random, itertools

random.seed(5)

def primes_upto(n):
    s = bytearray([1]) * (n + 1); s[0:2] = b'\x00\x00'
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]

def legendre(a, q):
    a %= q
    if a == 0:
        return 0
    return 1 if pow(a, (q - 1) // 2, q) == 1 else -1

def gp(p):
    for x in range(1, math.isqrt(p) + 1):
        y = math.isqrt(p - x * x)
        if y * y == p - x * x:
            return (x, y)

def mul(a, b):
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])

def conj(a):
    return (a[0], -a[1])

def E(n, Q):
    b, r = divmod(n, Q)
    return Q * b * (b - 1) // 2 + r * b

def gpow_mod(t, e, q):
    res = (1, 0); base = (t[0] % q, t[1] % q)
    while e:
        if e & 1:
            res = ((res[0] * base[0] - res[1] * base[1]) % q, (res[0] * base[1] + res[1] * base[0]) % q)
        base = ((base[0] * base[0] - base[1] * base[1]) % q, (2 * base[0] * base[1]) % q)
        e >>= 1
    return res

def run(ps, Xmax, nclusters, Msize):
    N = math.prod(ps)
    pis = [gp(p) for p in ps]
    units = [(1, 0), (0, 1), (-1, 0), (0, -1)]
    pts = []
    for a in itertools.product((0, 1), repeat=len(ps)):
        z = (1, 0)
        for j, aj in enumerate(a):
            z = mul(z, pis[j] if aj else conj(pis[j]))
        for e in range(4):
            pts.append((mul(units[e], z), a, e))
    assert all(z[0] ** 2 + z[1] ** 2 == N for z, _, _ in pts)
    pts.sort(key=lambda t: math.atan2(t[0][1], t[0][0]))
    ok = True; ncoll = 0; nchecks = 0
    clusters = [random.sample(pts, Msize) for _ in range(nclusters)]
    for s in range(0, len(pts) - Msize, max(1, len(pts) // nclusters)):
        clusters.append(pts[s:s + Msize])
    for q in primes_upto(Xmax):
        if q == 2 or N % q == 0:
            continue
        chi = 1 if q % 4 == 1 else -1
        nu = [1 if legendre(p, q) == -1 else 0 for p in ps]
        lam = 1 if legendre(2, q) == -1 else 0
        mq = q - chi
        invN = pow(N, -1, q)
        for cl in clusters:
            sig = [(sum(a[j] * nu[j] for j in range(len(ps))) + e * lam) % 2 for (_, a, e) in cl]
            # (1)
            for x in range(len(cl)):
                for y in range(x + 1, len(cl)):
                    t = mul(cl[x][0], conj(cl[y][0])); t = (t[0] * invN % q, t[1] * invN % q)
                    assert (t[0] ** 2 + t[1] ** 2) % q == 1
                    v = gpow_mod(t, mq // 2, q)
                    want = (1, 0) if sig[x] == sig[y] else (q - 1, 0)
                    ok &= (v == want); nchecks += 1
            # (2),(3) for each power
            Q, a = q, 1
            while Q <= Xmax:
                m = mq * q ** (a - 1)
                n1 = sum(sig); n0 = len(cl) - n1
                coll = 0
                for x in range(len(cl)):
                    for y in range(x + 1, len(cl)):
                        zx, zy = cl[x][0], cl[y][0]
                        if (zx[0] - zy[0]) % Q == 0 and (zx[1] - zy[1]) % Q == 0:
                            coll += 1
                            ok &= (sig[x] == sig[y])
                ok &= coll >= E(n0, m // 2) + E(n1, m // 2)
                ncoll += coll
                Q *= q; a += 1
    return ok, nchecks, ncoll

allok = True
for ps, Xmax, ncl, Ms in [([5, 13, 17, 29, 37], 60, 6, 40), ([5, 13, 17, 29, 37, 41], 80, 6, 60),
                          ([13, 17, 29, 37, 41, 53, 61], 90, 4, 80)]:
    ok, nc, ncoll = run(ps, Xmax, ncl, Ms)
    print(f"N = {'*'.join(map(str, ps))}: parity/Euler checks {nc}, collisions seen {ncoll}: {'OK' if ok else 'FAIL'}")
    allok &= ok
# (4)
ok4 = True
for Q in range(1, 13):
    for n0 in range(0, 40):
        for n1 in range(0, 40):
            ok4 &= E(n0, Q) + E(n1, Q) >= E(n0 + n1, 2 * Q)
            ok4 &= E(n0, Q) + E(n1, Q) <= E(n0 + n1, Q)
print("E-inequalities (4):", "OK" if ok4 else "FAIL")
print("ALL PARITY CHECKS PASSED" if (allok and ok4) else "SOME CHECK FAILED")
