# slverify.py -- INDEPENDENT checker for Rig^SL witnesses (round 15).  Imports nothing from base/ or the repository.
# Witness format (dict): P (primes), K, Q, F (forced set the witness must cover), phi ("genuine" or list of r floats
# given as exact decimal strings), and for every q in F: g1 = (a,b) a claimed generator of T_q (level 1), top = label
# vector mod Mtop = m_q q^(A-1), A = depth.
# Semantics (sharp.md Cor. 2.2 + Prop. 2.1): T_(q^a) is cyclic of order m_a = m_q q^(a-1), {1,i,-1,-i} is its unique
# subgroup of order 4, so  v collides at level a with unit s  <=>  m_a/4 | <v,ell>  and  <v,ell> = s*iota (mod m_q),
# iota = dlog_(g1)(i).  Requirement: no v in V_K collides at level A (then no lift collides deeper).
# Checks: (0) g1 generates T_q; (1) parity ell_j = [(p_j/q) = -1] mod 2; (2) coverage of F; (3) for every v,
# 1 <= ||v||_1 <= K (both signs) and s: m_(v,s) = prod q^(b_q) satisfies (1.2):
#     genuine phi:  m_(v,s) <= |D_s(Z_v)|                     (exact integers)
#     given phi :   nu_s m_(v,s) e^(-w/2) <= |sin(<v,phi> - s pi/4)|  in mpmath interval arithmetic (60 digits);
# (4) class status of every label: genuine class / some column-sign class / none (mod n_q, units of Z/n_q).
import sys, math, itertools, pickle, json
from mpmath import iv, mp

def legendre(a, q):
    t = pow(a % q, (q - 1) // 2, q)
    return -1 if t == q - 1 else t

def gp(p):
    for a in range(1, int(math.isqrt(p)) + 1, 2):
        b2 = p - a * a; b = math.isqrt(b2)
        if b * b == b2 and b % 2 == 0: return (a, b)
    raise ValueError

def gmul(x, y, q=None):
    z = (x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0])
    return (z[0] % q, z[1] % q) if q else z

def gpow(x, e, q):
    r = (1, 0)
    while e:
        if e & 1: r = gmul(r, x, q)
        x = gmul(x, x, q); e >>= 1
    return r

def tor(z, q):
    """z / conj z in Z[i]/q."""
    a, b = z[0] % q, z[1] % q; nn = (a * a + b * b) % q; inv = pow(nn, -1, q)
    return (((a * a - b * b) * inv) % q, (2 * a * b * inv) % q)

def factor_small(n):
    out = []; d = 2
    while d * d <= n:
        if n % d == 0:
            out.append(d)
            while n % d == 0: n //= d
        d += 1
    if n > 1: out.append(n)
    return out

def dlog_table(g, q, m):
    tab = {}; x = (1, 0)
    for k in range(m): tab[x] = k; x = gmul(x, g, q)
    return tab

def Zv(pis, v):
    z = (1, 0)
    for pi, t in zip(pis, v):
        f = pi if t > 0 else (pi[0], -pi[1])
        for _ in range(abs(t)): z = gmul(z, f)
    return z

def monos(r, K):
    return [v for v in itertools.product(range(-K, K + 1), repeat=r) if 0 < sum(map(abs, v)) <= K]

def check(W, verbose=True):
    P, K = W["P"], W["K"]; r = len(P); F = sorted(W["F"]); pis = [gp(p) for p in P]
    assert sorted(W["labels"]) == F, "witness must cover the forced set exactly"
    V = monos(r, K); mvs = {}; status = {}; epss = {}
    for q in F:
        lab = W["labels"][q]; g1 = tuple(lab["g1"]); top = list(lab["top"]); A = lab["A"]
        m = q - (1 if q % 4 == 1 else -1); n = m // 4; Mtop = m * q ** (A - 1)
        assert gmul(g1, (g1[0], -g1[1]), q) == (1, 0), "g1 not in T_q"
        assert all(gpow(g1, m // f, q) != (1, 0) for f in factor_small(m)), "g1 not a generator"
        tab = dlog_table(g1, q, m); iota = tab[(0, 1)]
        assert all(0 <= x < Mtop for x in top)
        for j, p in enumerate(P):
            assert top[j] % 2 == (1 if legendre(p, q) == -1 else 0), ("parity", q, j)
        for v in V:
            d = sum(a * b for a, b in zip(v, top)) % Mtop
            if d % n: continue
            s = [u for u in range(4) if (u * iota) % m == d % m]
            assert len(s) == 1; s = s[0]
            b = 1
            while b < A and d % (n * q ** b) == 0: b += 1
            assert not (b == A and d % (n * q ** (A - 1)) == 0 and A >= 1 and (d % (n * q ** (A - 1)) == 0) and b == A and A == b and False)
            if d % (n * q ** (A - 1)) == 0:
                raise AssertionError(f"collision at the deepest level A={A}, q={q}, v={v}")
            mvs[(v, s)] = mvs.get((v, s), 1) * q ** b
        L = [tab[tor(pi, q)] % n for pi in pis]; c = tuple(x % n for x in top)
        units = [u for u in range(1, n + 1) if math.gcd(u, n) == 1] if n > 1 else [1]
        E = set()
        for bits in range(2 ** r):
            eps = tuple(-1 if (bits >> j) & 1 else 1 for j in range(r))
            if any(tuple((u * e * x) % n for e, x in zip(eps, L)) == c for u in units): E.add(eps)
        epss[q] = E
        status[q] = "gen" if tuple([1] * r) in E or tuple([-1] * r) in E else ("eps" if E else "NONE")
    # (1.2)
    worst = None; nchk = 0
    if W["phi"] == "genuine":
        for v in V:
            z = Zv(pis, v); X, Y = z; D = (Y, X - Y, X, X + Y)
            for s in range(4):
                mm = mvs.get((v, s), 1); nchk += 1
                assert D[s] != 0 and mm <= abs(D[s]), ("slab at genuine phi", v, s, mm, D[s])
                sl = abs(D[s]) / mm
                worst = sl if worst is None or sl < worst else worst
    else:
        mp.dps = 60; iv.dps = 60
        phi = [iv.mpf(x) for x in W["phi"]]; pi4 = iv.pi / 4
        for v in V:
            N = 1
            for p, t in zip(P, v): N *= p ** abs(t)
            x = sum(t * f for t, f in zip(v, phi))
            for s in range(4):
                mm = mvs.get((v, s), 1); nchk += 1
                lhs = abs(iv.sin(x - s * pi4)); rhs = iv.mpf(mm) / iv.sqrt(iv.mpf(N)) * (iv.sqrt(2) if s % 2 else 1)
                gap = lhs - rhs
                assert gap.a > 0, ("slab at given phi", v, s, mm, float(lhs.a), float(rhs.b))
                g = float(gap.a); worst = g if worst is None or g < worst else worst
    # coherence
    seg = None; common = None
    for q in F:
        c2 = epss[q] if common is None else common & epss[q]
        if not c2: break
        common = c2; seg = q
    if verbose:
        print(f"OK: {len(F)} moduli, {nchk} (v,s) slab checks passed (phi={'genuine' if W['phi']=='genuine' else 'given'}); "
              f"min slack {'ratio |D|/m' if W['phi']=='genuine' else 'sin-gap'} = {worst}")
        print(f"   coherent (one eps) on forced q <= {seg};  first non-coherent: {[q for q in F if seg is None or q > seg][:1]}")
        print(f"   labels outside EVERY column-sign class: {sum(1 for q in F if status[q]=='NONE')}/{len(F)}; "
              f"in genuine class: {sum(1 for q in F if status[q]=='gen')}")
    return dict(seg=seg, status=status, worst=worst)

if __name__ == '__main__':
    W = pickle.load(open(sys.argv[1], "rb"))
    check(W)
