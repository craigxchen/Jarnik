"""Round 7, route G: tight monomial relations of genuine Gaussian primes and the pinning closure.

Exact integer arithmetic for every decision that enters a certificate.

Conventions (p3G.md section 2):
  * pi_j = a_j + i b_j, Norm = p_j (prime, = 1 mod 4); theta_j = arg pi_j.
  * A monomial v in Z^r has A_v = prod pi_j^(v_j^+) conj(pi_j)^(v_j^-) = X + iY.
  * D_s(A) = Y, X - Y, X, X + Y for s = 0, 1, 2, 3;  nu_s^2 = 1 (s even), 1/2 (s odd).
  * m_(v,s) = prod over odd q <= X, q not dividing N, of q^min(v_q(D_s), A_q),  A_q = max{a: q^a <= X}.
  * (v,s) is TIGHT iff m_(v,s) = |D_s(A_v)|, i.e. |D_s| is odd, coprime to N, and divides
    L_X = prod_{odd q <= X, q not | N} q^(A_q).
  * For a tight relation: S = sin^2 d = nu_s^2 D_s^2 / Norm(A_v)  (d = dist(arg A_v, s pi/4 + pi Z)),
    orientation sigma = sign sin(2(arg A_v - s pi/4)):
        s=0: sign(XY), s=1: sign(Y^2 - X^2), s=2: -sign(XY), s=3: sign(X^2 - Y^2).
    The admissible set for x = <v, phi - theta> is sigma*[0, pi - 2d] + pi Z.
Relations used: singles v = e_j, doubles v = 2 e_j, pairs v = e_i + tau e_j (i < j, tau = +-1).
"""
from fractions import Fraction
import numpy as np


def primes_upto(n):
    s = bytearray([1]) * (n + 1)
    s[0:2] = b'\x00\x00'
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]


def is_prime(n):
    if n < 2:
        return False
    small = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41)
    for p in small:
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2; s += 1
    for a in small:
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def isqrt(n):
    import math
    return math.isqrt(n)


def two_squares(p):
    """p prime = 1 mod 4 -> (a, b), a^2 + b^2 = p, a odd > 0, b even > 0 (Cornacchia/Hermite-Serret)."""
    # find x with x^2 = -1 mod p
    for c in range(2, p):
        if pow(c, (p - 1) // 2, p) == p - 1:
            x = pow(c, (p - 1) // 4, p)
            break
    a, b = p, x
    lim = isqrt(p)
    while b > lim:
        a, b = b, a % b
    u = b
    v = isqrt(p - u * u)
    assert u * u + v * v == p
    if u % 2 == 0:
        u, v = v, u
    return u, v


def lcm_odd_part(X, exclude=()):
    """L_X = prod over odd primes q <= X, q not in exclude, of q^(A_q)."""
    L = 1
    for q in primes_upto(X):
        if q == 2 or q in exclude:
            continue
        Q = q
        while Q * q <= X:
            Q *= q
        L *= Q
    return L


def smooth_mask(D, X, exclude=()):
    """Vectorised: True where |D| (int64 array, nonzero) divides L_X (exclude = primes of N <= X)."""
    g = np.abs(D.astype(np.int64))
    bad = np.zeros(g.shape, dtype=bool)
    for q in primes_upto(X):
        if q == 2:
            bad |= (g % 2 == 0)
            continue
        if q in exclude:
            bad |= (g % q == 0)
            continue
        A = 0
        Q = q
        while Q <= X:
            A += 1
            Q *= q
        for _ in range(A):
            m = (g % q == 0)
            g[m] //= q
        bad |= (g % q == 0)
    return (~bad) & (g == 1)


def sigma_of(Xv, Yv, s):
    """Orientation sign, exact (Python ints)."""
    def sg(t):
        return (t > 0) - (t < 0)
    Xv = int(Xv); Yv = int(Yv)
    if s == 0:
        return sg(Xv * Yv)
    if s == 1:
        return sg(Yv * Yv - Xv * Xv)
    if s == 2:
        return -sg(Xv * Yv)
    return sg(Xv * Xv - Yv * Yv)


def Ds(Xv, Yv, s):
    return (Yv, Xv - Yv, Xv, Xv + Yv)[s]


def tight_relations(gp, X, include_pairs=True, exclude=()):
    """gp: list of (a, b) Gaussian primes (the chosen pi_j). Returns list of relation dicts.

    Each relation: (kind, i, j, tau, s, D, normA, sigma)
      kind 'S' single v=e_i ; 'T' double v=2e_i ; 'P' pair v=e_i+tau e_j.
    All with exact Python ints for D, normA.  Uses int64 for the screening (|values| < 9e18 required).
    """
    a = np.array([g[0] for g in gp], dtype=np.int64)
    b = np.array([g[1] for g in gp], dtype=np.int64)
    norms = [g[0] * g[0] + g[1] * g[1] for g in gp]
    r = len(gp)
    rels = []
    # singles: A = pi_j
    for s in range(4):
        D = Ds(a, b, s)
        m = smooth_mask(D, X, exclude)
        for j in np.nonzero(m)[0]:
            j = int(j)
            Xv, Yv = int(a[j]), int(b[j])
            Dv = int(Ds(Xv, Yv, s))
            sg = int(sigma_of(Xv, Yv, s))
            rels.append(('S', j, j, 0, s, Dv, norms[j], sg))
    # doubles: A = pi_j^2 = (a^2 - b^2) + 2ab i
    X2 = a * a - b * b
    Y2 = 2 * a * b
    for s in range(4):
        D = Ds(X2, Y2, s)
        m = smooth_mask(D, X, exclude)
        for j in np.nonzero(m)[0]:
            j = int(j)
            Xv, Yv = int(X2[j]), int(Y2[j])
            Dv = int(Ds(Xv, Yv, s))
            sg = int(sigma_of(Xv, Yv, s))
            rels.append(('T', j, j, 0, s, Dv, norms[j] ** 2, sg))
    if not include_pairs:
        return rels
    iu, ju = np.triu_indices(r, k=1)
    ai, bi, aj, bj = a[iu], b[iu], a[ju], b[ju]
    for tau in (1, -1):
        # tau=+1: pi_i * pi_j ; tau=-1: pi_i * conj(pi_j)
        Xp = ai * aj - tau * bi * bj
        Yp = ai * (tau * bj) + bi * aj
        if tau == -1:
            Yp = -ai * bj + bi * aj
        for s in range(4):
            D = Ds(Xp, Yp, s)
            m = smooth_mask(D, X, exclude)
            for k in np.nonzero(m)[0]:
                i, j = int(iu[k]), int(ju[k])
                Xv, Yv = int(Xp[k]), int(Yp[k])
                Dv = int(Ds(Xv, Yv, s))
                sg = int(sigma_of(Xv, Yv, s))
                rels.append(('P', i, j, tau, s, Dv, norms[i] * norms[j], sg))
    return rels


def S_value(rel):
    """sin^2 d as an exact Fraction."""
    kind, i, j, tau, s, D, normA, sg = rel
    if s % 2 == 0:
        return Fraction(D * D, normA)
    return Fraction(D * D, 2 * normA)


def vec_of(rel):
    """The monomial vector as a dict {column: coeff}."""
    kind, i, j, tau, s, D, normA, sg = rel
    if kind == 'Q':
        (i0, j0, k0), (t1, t2) = i, j
        return {i0: 1, j0: t1, k0: t2}
    if kind == 'S':
        return {i: 1}
    if kind == 'T':
        return {i: 2}
    return {i: 1, j: tau}


# ------------------------------------------------------------------ pinning closure
class Pin:
    """Union-find encoding a sublattice Pi of Z^r known to satisfy <v, phi - theta> in pi Z.

    For every column k: e_k - sign[k] e_root(k) in Pi.  Root status: 0 free, 1 half (2 e_R in Pi),
    2 exact (e_R in Pi).
    """

    def __init__(self, r):
        self.par = list(range(r))
        self.sgn = [1] * r
        self.st = [0] * r

    def find(self, k):
        path = []
        while self.par[k] != k:
            path.append(k)
            k = self.par[k]
        root = k
        # path[-1] is a child of root; walk back accumulating sign relative to root
        acc = 1
        for node in reversed(path):
            acc *= self.sgn[node]
            self.sgn[node] = acc
            self.par[node] = root
        return root

    def rs(self, k):
        R = self.find(k)
        return R, (1 if R == k else self.sgn[k])

    def cls(self, vec):
        """Reduced class of vec (dict col->coeff) in Z^r / Pi: canonical tuple."""
        c = {}
        for k, co in vec.items():
            R, s = self.rs(k)
            c[R] = c.get(R, 0) + s * co
        out = []
        for R, co in c.items():
            st = self.st[R]
            if st == 2:
                continue
            if st == 1:
                co = co % 2
                if co:
                    out.append((R, 1))
                continue
            if co:
                out.append((R, co))
        return tuple(sorted(out))

    def neg(self, key):
        return tuple((R, co if self.st[R] == 1 else -co) for R, co in key)

    def add(self, key):
        """Add the class vector key (already reduced) to Pi. Returns True if Pi grew."""
        if not key:
            return False
        if len(key) == 1:
            (R, co), = key
            if abs(co) == 1:
                if self.st[R] < 2:
                    self.st[R] = 2
                    return True
                return False
            if abs(co) == 2:
                if self.st[R] < 1:
                    self.st[R] = 1
                    return True
                return False
            raise ValueError(key)
        if len(key) == 2:
            (R1, c1), (R2, c2) = key
            if abs(c1) != 1 or abs(c2) != 1:
                raise ValueError(key)
            # c1 e_R1 + c2 e_R2 in Pi  ->  e_R1 = -(c2/c1) e_R2 = eps e_R2
            eps = -c2 * c1
            self.par[R1] = R2
            self.sgn[R1] = eps
            self.st[R2] = max(self.st[R2], self.st[R1])
            return True
        raise ValueError(key)


def closure(r, rels, max_pass=200, verbose=False):
    """Pinning closure with circuits of length <= 2 (exact certificates).

    Returns (Pin object, log list of certificates).
    Certificate types:
      ('self', e)      : 2 sigma_e v_e in Pi and S_e > 1/2   (2 a_e < pi)
      ('pair', e, e')  : sigma_e v_e + sigma_e' v_e' in Pi and S_e + S_e' > 1  (a_e + a_e' < pi)
    """
    P = Pin(r)
    Sv = [S_value(rel) for rel in rels]
    vecs = [vec_of(rel) for rel in rels]
    alive = set(range(len(rels)))
    log = []
    for it in range(max_pass):
        grew = False
        buckets = {}
        for e in list(alive):
            key = P.cls(vecs[e])
            if not key:
                alive.discard(e)
                continue
            sg = rels[e][7]
            if sg == -1:
                key = P.neg(key)
            buckets.setdefault(key, []).append(e)
        done_keys = set()
        for key, lst in buckets.items():
            if key in done_keys:
                continue
            nk = P.neg(key)
            lst.sort(key=lambda e: Sv[e], reverse=True)
            cert = None
            if nk == key:
                if Sv[lst[0]] > Fraction(1, 2):
                    cert = ('self', lst[0])
                elif len(lst) >= 2 and Sv[lst[0]] + Sv[lst[1]] > 1:
                    cert = ('pair', lst[0], lst[1])
            elif nk in buckets:
                l2 = sorted(buckets[nk], key=lambda e: Sv[e], reverse=True)
                if Sv[lst[0]] + Sv[l2[0]] > 1:
                    cert = ('pair', lst[0], l2[0])
            done_keys.add(key); done_keys.add(nk)
            if cert is None:
                continue
            # recompute the class now (Pi may have grown in this pass)
            e = cert[1]
            k2 = P.cls(vecs[e])
            if not k2:
                continue
            # re-verify the certificate against the current Pi
            ok = False
            if cert[0] == 'self':
                kk = k2 if rels[e][7] == 1 else P.neg(k2)
                ok = (P.neg(kk) == kk)
            else:
                e2 = cert[2]
                ka = P.cls(vecs[e]); kb = P.cls(vecs[e2])
                ka = ka if rels[e][7] == 1 else P.neg(ka)
                kb = kb if rels[e2][7] == 1 else P.neg(kb)
                ok = (P.neg(ka) == kb)
            if not ok:
                continue
            try:
                added = P.add(k2)
            except ValueError:
                added = False  # class not representable in the union-find; skip (sound)
            if added:
                grew = True
                log.append(cert)
        if verbose:
            print('pass', it, 'alive', len(alive), 'certs', len(log))
        if not grew:
            break
    return P, log


def pinned_status(P, r):
    """Per column: 0 free, 1 pinned mod pi/2 (2 e_k in Pi), 2 pinned mod pi (e_k in Pi)."""
    out = []
    for k in range(r):
        R, s = P.rs(k)
        out.append(P.st[R])
    return out


def components(P, r):
    from collections import Counter
    return Counter(P.find(k) for k in range(r))


def triple_relations(gp, X, exclude=()):
    """Tight relations v = e_i + t1 e_j + t2 e_k (i<j<k, t1,t2 = +-1). int64 screening; needs |values| < 9e18."""
    a = np.array([g[0] for g in gp], dtype=np.int64)
    b = np.array([g[1] for g in gp], dtype=np.int64)
    norms = [g[0] * g[0] + g[1] * g[1] for g in gp]
    r = len(gp)
    rels = []
    for i in range(r):
        for j in range(i + 1, r):
            ks = np.arange(j + 1, r)
            if len(ks) == 0:
                continue
            for t1 in (1, -1):
                # A_ij = pi_i * pi_j^(t1)
                xi, yi = int(a[i]), int(b[i]); xj, yj = int(a[j]), t1 * int(b[j])
                X1, Y1 = xi * xj - yi * yj, xi * yj + yi * xj
                for t2 in (1, -1):
                    xk, yk = a[ks], t2 * b[ks]
                    Xp = X1 * xk - Y1 * yk
                    Yp = X1 * yk + Y1 * xk
                    for s in (1, 2, 3):
                        D = Ds(Xp, Yp, s)
                        m = smooth_mask(D, X, exclude)
                        for idx in np.nonzero(m)[0]:
                            k = int(ks[idx])
                            Xv, Yv = int(Xp[idx]), int(Yp[idx])
                            Dv = int(Ds(Xv, Yv, s))
                            sg = int(sigma_of(Xv, Yv, s))
                            rels.append(('Q', (i, j, k), (t1, t2), 0, s, Dv, norms[i] * norms[j] * norms[k], sg))
    return rels
