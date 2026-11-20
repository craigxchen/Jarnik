"""Referee check (independent code): factored witnesses r_j = i^{u_j} pi_j prod gamma_l^{V_jl},
local marks per sharp.md Def 1.1(d), exterior pair/character identities, same-level identities at
parent primes with gamma_l omitted, global collision divisibility.  Exact arithmetic in Z[i]/Q."""
import random, sys
from math import isqrt
from itertools import combinations
from sympy import primerange, factorint
sys.path.insert(0, '.')
from linalg_check import run, chi

def gmul(x, y, Q):
    return ((x[0]*y[0] - x[1]*y[1]) % Q, (x[0]*y[1] + x[1]*y[0]) % Q)
def gconj(x): return (x[0], -x[1])
def gnorm(x): return x[0]*x[0] + x[1]*x[1]
def ginv(x, Q):
    n = gnorm(x) % Q
    ni = pow(n, -1, Q)  # raises if not a unit
    return ((x[0]*ni) % Q, (-x[1]*ni) % Q)
def gpow(x, e, Q):
    x = (x[0] % Q, x[1] % Q)
    if e < 0: x, e = ginv(x, Q), -e
    r = (1 % Q, 0)
    while e:
        if e & 1: r = gmul(r, x, Q)
        x = gmul(x, x, Q); e >>= 1
    return r
IU = (0, 1)

def gauss_factor(p):
    for x in range(1, isqrt(p)+1):
        y2 = p - x*x; y = isqrt(y2)
        if y*y == y2 and x > y: return (x, y)   # argument in (0, pi/4)
    raise ValueError

def check(q, rho, seed, winding_mode="even", ext_primes=None, powers=(1, 2, 3), nchar=200, verbose=True):
    data = run(q, rho, seed, verbose=False, do_snf=False)
    labels, S, a, U, V, F, i0 = data['labels'], data['S'], data['a'], data['U'], data['V'], data['F'], data['i0']
    n = len(labels); m = len(S)
    rng = random.Random(1000 + seed)
    primes = [p for p in primerange(5, 10**5) if p % 4 == 1][:n]
    rng.shuffle(primes)                       # random assignment of primes to columns
    pis = [gauss_factor(p) for p in primes]   # pi_j for label labels[j]
    col = {lab: jj for jj, lab in enumerate(labels)}
    # windings
    if winding_mode == "even":
        k = {i: 2*rng.randrange(-5, 6) + 1 for i in S}       # all odd, differences even, nonconstant
    else:
        k = {i: rng.randrange(-5, 6) for i in S}
        if all((k[i]-k[i0]) % 2 == 0 for i in S): k[S[-1]] += 1
    others = [i for i in S if i != i0]
    if winding_mode == "even":
        rhs = [(k[i] - k[i0]) // 2 for i in others]
    else:
        rhs = [(k[i] - k[i0]) // 2 for i in others]  # floor: deliberately wrong bookkeeping
    u = [sum(U[jj][c] * rhs[c] for c in range(m-1)) for jj in range(n)]
    if ext_primes is None:
        ext_primes = [p for p in primerange(3, 80) if p not in primes]

    def witnesses(Q, omit=None):
        gam = [None if jj == omit else gmul(pis[jj], ginv(gconj(pis[jj]), Q), Q) for jj in range(n)]
        r = []
        for jj in range(n):
            if jj == omit: r.append(None); continue
            x = gmul(gpow(IU, u[jj], Q), pis[jj], Q)
            for ll in range(n):
                if ll != omit and V[jj][ll]:
                    x = gmul(x, gpow(gam[ll], V[jj][ll], Q), Q)
            assert gnorm(x) % Q == primes[jj] % Q          # (d): Norm r_j = p_j
            r.append(x)
        return gam, r

    def marks(Q, r, omit=None):
        rho_ = {}
        for i in S:
            x = gpow(IU, -k[i], Q)
            for jj in range(n):
                if jj == omit: continue
                e = a[i][jj]
                x = gmul(x, gmul(gpow(r[jj], e, Q), gpow(gconj(r[jj]), 2 - e, Q), Q), Q)
            rho_[i] = x
        return rho_

    ext_pairs = ext_chars = 0; fails = 0
    for ell in ext_primes:
        for pw in powers:
            Q = ell**pw
            gam, r = witnesses(Q)
            # check r_j/conj r_j = (-1)^{u_j} prod gamma^{F_jl}
            for jj in range(n):
                lhs = gmul(r[jj], ginv(gconj(r[jj]), Q), Q)
                rr = gpow((Q-1, 0), u[jj] % 2, Q)
                for ll in range(n):
                    if F[jj][ll]: rr = gmul(rr, gpow(gam[ll], F[jj][ll], Q), Q)
                assert lhs == rr
            rh = marks(Q, r)
            for i, i2 in combinations(S, 2):
                ok = gmul(rh[i], gam[col[i2]], Q) == gmul(rh[i2], gam[col[i]], Q)
                ext_pairs += 1; fails += (not ok)
            for _ in range(nchar):
                c = {i: rng.randrange(-3, 4) for i in S}
                c[i0] -= sum(c.values())
                lhs = (1, 0); rhs_ = (1, 0)
                for i in S:
                    lhs = gmul(lhs, gpow(rh[i], c[i], Q), Q)
                    rhs_ = gmul(rhs_, gpow(gam[col[i]], c[i], Q), Q)
                ext_chars += 1; fails += (lhs != rhs_)
    # parent primes
    lev_pairs = 0; diff_level_match = diff_level_total = 0; owner_unique = True
    for ll, lab in enumerate(labels):
        p = primes[ll]
        levels = {i: a[i][ll] for i in S}
        # owner: unique mark with level 1
        ones = [i for i in S if levels[i] == 1]
        if lab in S:
            owner_unique &= (ones == [lab])
        else:
            owner_unique &= (ones == [])
        for pw in powers:
            Q = p**pw
            gam, r = witnesses(Q, omit=ll)
            rh = marks(Q, r, omit=ll)
            for i, i2 in combinations(S, 2):
                if levels[i] == levels[i2]:
                    assert lab not in (i, i2)
                    ok = gmul(rh[i], gam[col[i2]], Q) == gmul(rh[i2], gam[col[i]], Q)
                    lev_pairs += 1; fails += (not ok)
                elif lab not in (i, i2):
                    diff_level_total += 1
                    diff_level_match += gmul(rh[i], gam[col[i2]], Q) == gmul(rh[i2], gam[col[i]], Q)
    if verbose:
        print(f"q={q} rho={rho} seed={seed} m={m} windings={winding_mode}: exterior pair identities={ext_pairs}, "
              f"exterior random characters={ext_chars}, same-level identities={lev_pairs}, failures={fails}; "
              f"owner unique={owner_unique}; different-level pairs with the same ratio={diff_level_match}/{diff_level_total}")
    return fails, dict(primes=primes, pis=pis, S=S, a=a, labels=labels, col=col, k=k, u=u, V=V)

def global_collision(q, rho, seed, bound=200):
    """For pairs: collision of gamma_i/gamma_i' with unit u0 at Q (odd, prime to p_i p_i') iff Q | A - u0 conj A,
    A = pi_i conj(pi_i'); and D_1 <= L0 |beta_i - beta_i'|."""
    import cmath, math
    _, d = check(q, rho, seed, verbose=False, nchar=0, ext_primes=[3])
    pis, S, col, primes = d['pis'], d['S'], d['col'], d['primes']
    units = [(1, 0), (0, 1), (-1, 0), (0, -1)]
    checked = 0; maxratio = 0
    L0 = max(primes[col[i]] for i in S)
    for i, i2 in combinations(S, 2):
        x, y = pis[col[i]], pis[col[i2]]
        A = (x[0]*y[0] + x[1]*y[1], x[1]*y[0] - x[0]*y[1])     # x * conj(y)
        for u0 in units:
            uc = (u0[0]*A[0] + u0[1]*A[1], u0[1]*A[0] - u0[0]*A[1])   # u0 * conj(A)
            Dg = (A[0] - uc[0], A[1] - uc[1])
            assert Dg != (0, 0)
            assert (Dg[0] + Dg[1]) % 2 == 0                          # divisible by 1+i
            for ell in primerange(3, bound):
                if ell in (primes[col[i]], primes[col[i2]]): continue
                for pw in (1, 2):
                    Q = ell**pw
                    g1 = gmul(x, ginv(gconj(x), Q), Q); g2 = gmul(y, ginv(gconj(y), Q), Q)
                    coll = gmul(g1, ginv(g2, Q), Q) == (u0[0] % Q, u0[1] % Q)
                    div = Dg[0] % Q == 0 and Dg[1] % Q == 0
                    assert coll == div; checked += 1
        Im = A[1]
        b1 = math.atan2(x[1], x[0]); b2 = math.atan2(y[1], y[0])
        odd = abs(Im)
        while odd % 2 == 0: odd //= 2
        assert odd <= L0 * abs(b1 - b2) + 1e-9
        assert abs(b1 - b2) >= 1 / L0
        maxratio = max(maxratio, odd / (L0 * abs(b1 - b2)))
    print(f"global collision test q={q}: {checked} (pair, unit, Q) equivalences collision <=> Q | A - u0 conj A; "
          f"max odd|Im|/(L0|dbeta|) = {maxratio:.4f}")

if __name__ == "__main__":
    tot = 0
    for q, rho, seed in [(11, 1.0, 0), (19, 1.0, 0), (19, 0.8, 1), (23, 0.9, 2), (31, 0.9, 0)]:
        f, _ = check(q, rho, seed, powers=(1, 2, 3) if q <= 19 else (1, 2))
        tot += f
    f, _ = check(19, 1.0, 0, winding_mode="odd", nchar=20)
    print("negative control (an odd winding difference) failures:", f)
    global_collision(19, 1.0, 0)
    print("TOTAL failures in positive tests:", tot)
