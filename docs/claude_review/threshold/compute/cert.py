"""Exact lower-bound certificates  reg_lam(A) >= r  (hence reg(A) >= r).

Given (M, point set A = shell or slice, lam, r):
 1. kernel of the generator matrix of degree <= r-1, mod three primes (canonical RREF kernel
    vector of the first free column), CRT + rational reconstruction -> integer gamma != 0;
 2. the measure x = a_t b_t sum_T gamma_T 1_{k_T} is nonzero because the c_T = a_t b_t 1_{k_T}
    are linearly independent (rank of the full generator matrix = N mod P1, which implies
    rank N over Q);
 3. x is orthogonal to ALL polynomials of degree <= r-1: since <a_t b_t z, P> = <z, b_t a_t P>
    it suffices that sum_T gamma_T (b_t a_t x^beta)(k_T) = 0 for every monomial x^beta of
    degree <= r-1 (beta up to the row group R_t, which leaves a_t x^beta unchanged).  These
    integers are bounded by B = sum|gamma| * |R||C| * t^(r-1); they are checked to vanish
    modulo primes whose product exceeds 2B, so they vanish over Z.
Step 3 uses only elementary facts (monomials span Poly_{<=r-1}); it does NOT use the
Ariki-Terasoma-Yamada spanning theorem.  Optionally (--direct) x is also built explicitly and
its moments in k_2..k_M are computed with Python integers.
"""
import sys, itertools, math
import numpy as np
from fractions import Fraction
from modla import Echelon, P1, P2, is_prime
from combi import partitions, shell_orbits, slice_orbits, hook_dim
from iso import LamData, columns_for, sym_filtration, subset_dp, _powers, group_perms

PRIMES = [P1, P2]
q = P2 - 2
while len(PRIMES) < 40:
    if is_prime(q):
        PRIMES.append(q)
    q -= 2


def F_general(ld, K, alphas, p):
    """(b_t a_t x^alpha)(k) mod p for each exponent vector alpha (rows) and point k (cols)."""
    saved = ld.alpha
    ld.alpha = [list(a) for a in alphas]
    ld.f_saved = ld.f
    ld.f = len(alphas)
    try:
        out = ld.F_values(K, p)
    finally:
        ld.alpha = saved
        ld.f = ld.f_saved
    return out


def gen_matrix(ld, orbits, p, dmax):
    """generator rows of degree <= dmax (as in iso.lam_profile), columns = SSYT points."""
    M = ld.M
    K, oidx, relevant = columns_for(ld.lam, orbits)
    V = np.array([orbits[o] for o in relevant], dtype=np.int64)
    remap = {o: i for i, o in enumerate(relevant)}
    oloc = np.array([remap[o] for o in oidx], dtype=np.int64)
    basis, _ = sym_filtration(V, M, p)
    G = np.array([b[1] for b in basis], dtype=np.int64)
    ge = np.array([b[0] for b in basis], dtype=np.int64)
    Fv = ld.F_values(K, p)
    rows = []
    for s, c in enumerate(ld.charge):
        for j in np.flatnonzero(ge + c <= dmax):
            rows.append((G[j][oloc] * Fv[s]) % p)
    B = np.array(rows, dtype=np.int64).reshape(-1, K.shape[0])
    return B, K


def ratrec(a, m):
    bound = math.isqrt(m // 2)
    r0, r1 = m, a % m
    s0, s1 = 0, 1
    while r1 > bound:
        qq = r0 // r1
        r0, r1 = r1, r0 - qq * r1
        s0, s1 = s1, s0 - qq * s1
    if s1 == 0 or abs(s1) > bound:
        return None
    return Fraction(r1, s1)


def kernel_vector(ld, orbits, r, nprimes=3, maxprimes=12):
    """integer kernel vector gamma of the degree <= r-1 generator matrix (or None)."""
    images = []
    for p in PRIMES[:maxprimes]:
        B, K = gen_matrix(ld, orbits, p, r - 1)
        E = Echelon(K.shape[0], p)
        if B.shape[0]:
            E.add(B)
        Kb = E.kernel_basis()
        if Kb.shape[0] == 0:
            return None, K, 'no kernel mod %d' % p
        images.append((p, tuple(sorted(int(x) for x in E.piv)), Kb[0]))
        if len(images) < nprimes:
            continue
        pivs = set(x[1] for x in images)
        if len(pivs) != 1:
            return None, K, 'pivot structure differs between primes'
        Mod = 1
        for (pp, _, _) in images:
            Mod *= pp
        n = K.shape[0]
        vals = []
        ok = True
        for j in range(n):
            res = 0
            for (pp, _, v) in images:
                Mp = Mod // pp
                res = (res + int(v[j]) * Mp * pow(Mp, -1, pp)) % Mod
            fr = ratrec(res, Mod)
            if fr is None:
                ok = False
                break
            vals.append(fr)
        if not ok:
            continue          # add another prime
        L = 1
        for fr in vals:
            L = L * fr.denominator // math.gcd(L, fr.denominator)
        g = [int(fr * L) for fr in vals]
        gg = 0
        for x in g:
            gg = math.gcd(gg, x)
        g = [x // gg for x in g]
        return g, K, 'ok (%d primes)' % len(images)
    return None, None, 'rational reconstruction failed'


def beta_reps(ld, dmax):
    """representatives of R_t-orbits of exponent vectors beta in N^M with |beta| <= dmax."""
    M = ld.M
    rows = ld.rows
    out = []
    for d in range(dmax + 1):
        # distribute d among rows, then a partition (sorted tuple) within each row
        def rec(ri, rem, cur):
            if ri == len(rows):
                if rem == 0:
                    beta = [0] * M
                    for rr, part in zip(rows, cur):
                        for posn, e in zip(rr, part):
                            beta[posn] = e
                    out.append(beta)
                return
            L = len(rows[ri])
            for take in range(rem + 1):
                for part in _parts_len(take, L):
                    rec(ri + 1, rem - take, cur + [part])
        rec(0, d, [])
    return out


def _parts_len(n, L):
    """weakly decreasing L-tuples of nonnegative integers with sum n."""
    def rec(n, L, mx):
        if L == 0:
            if n == 0:
                yield ()
            return
        for a in range(min(n, mx), -1, -1):
            if a * L < n:
                break
            for rest in rec(n - a, L - 1, a):
                yield (a,) + rest
    return list(rec(n, L, n))


def verify(ld, K, gamma, r, t, chunk=4000):
    """check sum_T gamma_T (b a x^beta)(k_T) = 0 for all beta, |beta| <= r-1 (mod enough primes)."""
    betas = beta_reps(ld, r - 1)
    B = sum(abs(g) for g in gamma) * ld.nR * ld.nC * max(1, t) ** (r - 1)
    need = 2 * B
    prod = 1
    used = 0
    for p in PRIMES:
        gp = np.array([g % p for g in gamma], dtype=np.int64)
        for c0 in range(0, len(betas), chunk):
            Fv = F_general(ld, K, betas[c0:c0 + chunk], p)
            # dot products mod p without overflow: split gamma
            g1, g0 = gp >> 16, gp & 0xFFFF
            v = ((Fv * g1[None, :]) % p).sum(axis=1) % p
            v = (v * 65536 + ((Fv * g0[None, :]) % p).sum(axis=1)) % p
            if np.any(v != 0):
                return False, len(betas), used
        prod *= p
        used += 1
        if prod > need:
            return True, len(betas), used
    return None, len(betas), used


def direct_moments(ld, K, gamma, r, s):
    """build x = a_t b_t sum gamma_T 1_{k_T} explicitly and check all moments of degree <= r-1
    in k_2..k_M with Python integers (small cases only).  Returns (ok, support size)."""
    M = ld.M
    Cp, Cs = group_perms(ld.cols, M, True)
    Rp, _ = group_perms(ld.rows, M, False)
    x = {}
    for gT, k in zip(gamma, K.tolist()):
        if gT == 0:
            continue
        for perm, sg in zip(Cp.tolist(), Cs.tolist()):
            ks = [k[perm[i]] for i in range(M)]
            for rp in Rp.tolist():
                kk = tuple(ks[rp[i]] for i in range(M))
                x[kk] = x.get(kk, 0) + sg * gT
    x = {k: v for k, v in x.items() if v}
    if not x:
        return False, 0
    pts = list(x.keys())
    for d in range(r):
        for mon in itertools.combinations_with_replacement(range(1, M), d):
            tot = 0
            for k in pts:
                v = x[k]
                for i in mon:
                    v *= k[i]
                tot += v
            if tot != 0:
                return False, len(pts)
    return True, len(pts)


def certify(M, pp, nn, lam, r, kind='shell', direct=False):
    orbits = shell_orbits(M, pp, nn) if kind == 'shell' else slice_orbits(M, pp, nn)
    ld = LamData(lam)
    gamma, K, msg = kernel_vector(ld, orbits, r)
    if gamma is None:
        return {'ok': False, 'msg': msg}
    ok, nb, used = verify(ld, K, gamma, r, pp + nn)
    res = {'ok': ok, 'msg': msg, 'nbeta': nb, 'primes_used': used, 'N': K.shape[0],
           'max_abs_gamma': max(abs(g) for g in gamma), 'nnz': sum(1 for g in gamma if g)}
    if direct:
        dok, supp = direct_moments(ld, K, gamma, r, pp - nn)
        res['direct_ok'] = dok
        res['support'] = supp
    return res


if __name__ == '__main__':
    M, pp, nn, r = map(int, sys.argv[1:5])
    lam = tuple(map(int, sys.argv[5].split(',')))
    kind = sys.argv[6] if len(sys.argv) > 6 else 'shell'
    direct = len(sys.argv) > 7 and sys.argv[7] == 'direct'
    print(certify(M, pp, nn, lam, r, kind, direct))


def verify_gen(ld, orbits, gamma, r, t):
    """ATY-level check: the degree <= r-1 generator matrix annihilates gamma over Z (checked
    modulo primes whose product exceeds twice the bound on the entries of G gamma)."""
    M = ld.M
    B = sum(abs(g) for g in gamma) * ld.nR * ld.nC * (M * max(1, t)) ** (r - 1)
    prod, used = 1, 0
    for p in PRIMES:
        G, K = gen_matrix(ld, orbits, p, r - 1)
        gp = np.array([g % p for g in gamma], dtype=np.int64)
        g1, g0 = gp >> 16, gp & 0xFFFF
        v = ((G * g1[None, :]) % p).sum(axis=1) % p
        v = (v * 65536 + ((G * g0[None, :]) % p).sum(axis=1)) % p
        if np.any(v != 0):
            return False, used
        prod *= p
        used += 1
        if prod > 2 * B:
            return True, used
    return None, used


def count_betas(ld, dmax):
    from math import comb
    # rough count of R-orbit representatives
    tot = comb(dmax + ld.M, ld.M)
    return tot // max(1, ld.nR)


def certify2(M, pp, nn, lam, r, kind='shell', beta_limit=300000):
    orbits = shell_orbits(M, pp, nn) if kind == 'shell' else slice_orbits(M, pp, nn)
    ld = LamData(lam)
    if r <= 0:
        return {'ok': True, 'level': 'trivial', 'msg': 'r <= 0'}
    gamma, K, msg = kernel_vector(ld, orbits, r)
    if gamma is None:
        return {'ok': False, 'msg': msg}
    res = {'msg': msg, 'N': K.shape[0], 'max_abs_gamma': max(abs(g) for g in gamma),
           'nnz': sum(1 for g in gamma if g)}
    if count_betas(ld, r - 1) <= beta_limit:
        ok, nb, used = verify(ld, K, gamma, r, pp + nn)
        res.update(ok=ok, level='elementary (all monomials)', nbeta=nb, primes_used=used)
    else:
        ok, used = verify_gen(ld, orbits, gamma, r, pp + nn)
        res.update(ok=ok, level='ATY generators', primes_used=used)
    return res
