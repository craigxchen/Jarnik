"""EXACT (multimodular) isotypic computation of Psi_A on the Fourier circle C_1 of the sdc_q slab box A = A_L.

For m < q:  Psi_A(u_theta) = 4^{-m} sum_{i=0}^m binom(m,i)^2 ||R_i||^2  (independent of theta), where
F_i = khat_1^i khat_{-1}^{m-i} (component c_i = 2i - m of the translation group), R_i its L^2(mu_A)-residual modulo
Poly_{<m}|_A (which lies in the same component), and ||R_i||^2 is RATIONAL (AGL(1,q)/Galois invariance).
Also E_{mu_A} l_{u_theta}^{2m} = 4^{-m} sum binom(m,i)^2 ||F_i||^2 (constant), and the per-component dimensions
dim Poly_{<m,c} (exact mod p; rank agreement over all primes used is recorded).

Arithmetic: for primes p = 1 mod q (zeta -> z in F_p), complex conjugation = (z -> z^{-1}).  Column set of
Poly_{<m,c}: P_c(k^alpha) over monomial orbit representatives alpha (exponents <= 2) of degree < m, restricted to
monomials whose support fits in the max support of A (others vanish on A).  ||R||^2 = (1/U)(F*WF - F*WM G^{-1} M*WF),
G = M*WM on an independent column subset S (pivots mod p; G invertible mod p certifies independence over Q(zeta)).
Rational reconstruction from CRT over primes; reported once stable over 3 further primes.

usage: python3 exact_iso.py q L mmax [nprimes_max]
"""
import sys, time, itertools, math
import numpy as np
import flint
from fractions import Fraction

q, L, mmax = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
NPMAX = int(sys.argv[4]) if len(sys.argv) > 4 else 400
t0 = time.time()
reps = np.fromfile(f"data/A_q{q}_L{L}.bin", dtype=np.int8).reshape(-1, q)
fixed = [r for r in reps if (r == r[0]).all()]
free = np.array([r for r in reps if not (r == r[0]).all()], dtype=np.int64)
assert len(fixed) == 1 and not fixed[0].any()
assert len(reps) * (2 ** 24) ** 2 < 2 ** 63 and q * (2 ** 24) ** 2 < 2 ** 63, 'int64 overflow guard (p < 2^24)'
R = len(free); U = q * R + 1
smaxA = int((free != 0).sum(axis=1).max())
print(f"sdc{q} L={L}: |A|={U} free orbits {R}, max support {smaxA}", flush=True)
# all shifts of reps: orb[s] = g^s rep, (g^s k)_t = k_{t-s}
orb = np.stack([np.roll(free, s, axis=1) for s in range(q)])        # q x R x q


def mono_reps(deg):
    """orbit representatives (under cyclic shift) of exponent vectors in {0,1,2}^q, |a| = deg, support <= smaxA."""
    out = []; seen = set()
    for supp_size in range(0, min(deg, smaxA) + 1):
        if supp_size * 2 < deg: continue
        for S in itertools.combinations(range(q), supp_size):
            if supp_size and S[0] != 0: continue          # rotate so that 0 is in the support
            # exponents on S: each 1 or 2, sum deg
            for ex in itertools.product((1, 2), repeat=supp_size):
                if sum(ex) != deg: continue
                a = [0] * q
                for x, e in zip(S, ex): a[x] = e
                a = tuple(a)
                if a in seen: continue
                for s in range(q): seen.add(a[-s:] + a[:-s] if s else a)
                out.append(a)
    return out


MONO = {}
for e in range(mmax):
    MONO[e] = mono_reps(e)
print("monomial orbit reps by degree:", [len(MONO[e]) for e in range(mmax)], f"({time.time()-t0:.0f}s)", flush=True)


def mono_vals(alist):
    """integer array (q, R, len) of monomial values at all shifts of the free reps."""
    A = np.array(alist, dtype=np.int64).reshape(-1, q)        # N x q
    out = np.ones((q, R, len(alist)), dtype=np.int64)
    for x in range(q):
        ex = A[:, x]
        if not ex.any(): continue
        v = orb[:, :, x][:, :, None]                          # q x R x 1
        out *= np.where(ex[None, None, :] == 0, 1, np.where(ex[None, None, :] == 1, v, v * v))
    return out


MV = {e: mono_vals(MONO[e]) for e in range(mmax)}
print(f"monomial values computed ({time.time()-t0:.0f}s)", flush=True)
# fixed point k = 0: monomial value 1 iff degree 0


def primes_1modq(start, count):
    out = []; p = start - (start % q) + 1
    while len(out) < count:
        if p > start and flint.fmpz(p).is_prime(): out.append(p)
        p += q
    return out


PIV = {}


def run_prime(p, m, full_rank_check=False):
    """returns dict c -> rank of Poly_{<m,c}, and list over i of (||R_i||^2 * U mod p, ||F_i||^2 * U mod p)."""
    g0 = 2
    while pow(g0, (p - 1) // q, p) == 1: g0 += 1
    z = pow(g0, (p - 1) // q, p); zi = pow(z, p - 2, p)
    W = q % p
    # khat_1, khat_{-1} at reps
    zp = np.array([pow(z, t, p) for t in range(q)], dtype=np.int64)
    zm = np.array([pow(zi, t, p) for t in range(q)], dtype=np.int64)
    kh1 = (free % p) @ zp % p; khm = (free % p) @ zm % p
    res = {}; ranks = {}
    for i in range(m + 1):
        c = (2 * i - m) % q
        # columns: P_c(k^alpha) = sum_s z^{-cs} (g^s k)^alpha ; conj columns with z -> z^{-1}
        wz = np.array([pow(zi, (c * s) % q, p) for s in range(q)], dtype=np.int64)
        wzb = np.array([pow(z, (c * s) % q, p) for s in range(q)], dtype=np.int64)
        cols = []; colsb = []
        for e in range(m):
            mv = MV[e] % p                                   # q x R x N
            cols.append(np.tensordot(wz, mv, axes=(0, 0)) % p)
            colsb.append(np.tensordot(wzb, mv, axes=(0, 0)) % p)
        Mx = np.concatenate(cols, axis=1); Mb = np.concatenate(colsb, axis=1)
        rows_extra = (c == 0)
        if rows_extra:   # fixed point 0: P_0(k^alpha)(0) = q * [deg alpha == 0]; weight 1
            fx = np.zeros((1, Mx.shape[1]), dtype=np.int64); fx[0, 0] = q % p
            Mx = np.vstack([Mx, fx]); Mb = np.vstack([Mb, fx])
        wvec = np.full(Mx.shape[0], W, dtype=np.int64)
        if rows_extra: wvec[-1] = 1
        F = (np.array([pow(int(a), i, p) * pow(int(b), m - i, p) % p for a, b in zip(kh1, khm)], dtype=np.int64))
        Fb = (np.array([pow(int(b), i, p) * pow(int(a), m - i, p) % p for a, b in zip(kh1, khm)], dtype=np.int64))
        if rows_extra:
            F = np.append(F, 0 if m > 0 else 1); Fb = np.append(Fb, 0 if m > 0 else 1)
        nr, nc = Mx.shape
        if (m, c) in PIV and not full_rank_check:
            piv = PIV[(m, c)]; rk = len(piv)
        else:
            Mt = flint.nmod_mat(nr, nc, [int(v) for v in Mx.ravel()], p)
            Rr, rk = Mt.rref()                      # pivot columns of the rref = an independent spanning subset
            piv = []; r_ = 0
            for col in range(nc):
                if r_ < rk and int(Rr[r_, col]) != 0:
                    piv.append(col); r_ += 1
            PIV.setdefault((m, c), piv)
        ranks[c] = rk
        FWF = int((Fb * wvec % p * F % p).sum() % p)
        if rk == 0:
            res[i] = (FWF, FWF); continue
        MS = Mx[:, piv]; MSb = Mb[:, piv]
        G = (MSb.T * wvec) % p
        b = G @ F % p                                            # M*WF  (rk)
        bb = (MS.T * wvec) % p @ Fb % p                          # (MS^T W Fb) = conj-side
        Gm = flint.nmod_mat(rk, rk, [int(v) for v in (G @ MS % p).ravel()], p)
        if Gm.det() == 0:
            raise RuntimeError(f"Gram singular mod {p} at c={c}")
        sol = Gm.solve(flint.nmod_mat(rk, 1, [int(v) for v in b], p))
        proj = sum(int(bb[t]) * int(sol[t, 0]) for t in range(rk)) % p
        res[i] = ((FWF - proj) % p, FWF)
    return ranks, res


def total_rank(p, m):
    """sum over all components c of rank Poly_{<m,c} mod p (plus nothing else): = dim Poly_{<m}|_A mod p."""
    if m <= 0: return 0
    g0 = 2
    while pow(g0, (p - 1) // q, p) == 1: g0 += 1
    z = pow(g0, (p - 1) // q, p); zi = pow(z, p - 2, p)
    tot = 0
    for c in range(q):
        wz = np.array([pow(zi, (c * s) % q, p) for s in range(q)], dtype=np.int64)
        Mx = np.concatenate([np.tensordot(wz, MV[e] % p, axes=(0, 0)) % p for e in range(m)], axis=1)
        if c == 0:
            fx = np.zeros((1, Mx.shape[1]), dtype=np.int64); fx[0, 0] = q % p; Mx = np.vstack([Mx, fx])
        tot += flint.nmod_mat(Mx.shape[0], Mx.shape[1], [int(v) for v in Mx.ravel()], p).rank()
    return tot


def ratrec(a, mod):
    """rational reconstruction of a mod `mod` with |num|,den <= sqrt(mod/2)."""
    a %= mod; bound = math.isqrt(mod // 2)
    r0, r1 = mod, a; s0, s1 = 0, 1
    while r1 > bound:
        qq = r0 // r1
        r0, r1 = r1, r0 - qq * r1; s0, s1 = s1, s0 - qq * s1
    if s1 == 0 or abs(s1) > bound: return None
    fr = Fraction(r1, s1)
    return fr


for m in range(1, mmax + 1):
    if m >= q: break
    t1 = time.time()
    primes = primes_1modq(2 ** 23, NPMAX)   # p < 2^24: (R+1)*p^2 < 2^63 checked above
    acc = {}; mod = 1; stable = {}; last = {}; ranks0 = None; done = False
    for npr, p in enumerate(primes):
        ranks, res = run_prime(p, m, full_rank_check=(npr < 2))
        if ranks0 is None: ranks0 = ranks
        elif ranks != ranks0:
            print(f"  rank mismatch at prime {p}: {ranks} vs {ranks0}", flush=True)
        # CRT update
        if npr == 0:
            crt = {k: v for k, v in [(('R', i), r[0]) for i, r in res.items()] + [(('F', i), r[1]) for i, r in res.items()]}
            mod = p
        else:
            newv = {('R', i): r[0] for i, r in res.items()}; newv.update({('F', i): r[1] for i, r in res.items()})
            for k in crt:
                a, b_ = crt[k], newv[k]
                # x = a mod mod, x = b_ mod p
                t = ((b_ - a) * pow(mod, -1, p)) % p
                crt[k] = a + mod * t
            mod *= p
        rec = {k: ratrec(v, mod) for k, v in crt.items()}
        for k in rec:
            if rec[k] is not None and last.get(k) == rec[k]: stable[k] = stable.get(k, 0) + 1
            else: stable[k] = 0
        last = rec
        if all(stable.get(k, 0) >= 3 for k in crt):
            done = True; break
    deficits = {c: (R + (1 if c == 0 else 0)) - ranks0[c] for c in sorted(ranks0)}
    totrank = total_rank(primes[0], m); dimV = U - totrank
    totrank1 = total_rank(primes[0], m - 1) if m >= 1 else 0
    hA = totrank - totrank1   # = dim Poly_{<m} - dim Poly_{<m-1} = h_A(m-1)
    if not done:
        print(f"m={m}: NOT stabilised after {npr+1} primes", flush=True); continue
    RR = {i: last[('R', i)] / U for i in range(m + 1)}       # ||R_i||^2 in mu_A normalisation
    FF = {i: last[('F', i)] / U for i in range(m + 1)}
    Psi = sum(Fraction(math.comb(m, i) ** 2) * RR[i] for i in range(m + 1)) / 4 ** m
    E2m = sum(Fraction(math.comb(m, i) ** 2) * FF[i] for i in range(m + 1)) / 4 ** m
    print(f"m={m}: primes used {npr+1}; dim Poly_<m|_A = {totrank}, dim V(A,m) = {dimV}, h_A(m-1) = {hA}; deficits delta_c on needed components {deficits}", flush=True)
    for i in range(m + 1):
        print(f"   i={i} c={(2*i-m)%q}: ||R_i||^2 = {RR[i]}  (~{float(RR[i]):.6e});  ||F_i||^2 = {FF[i]} (~{float(FF[i]):.6e})")
    ratio = Psi / E2m if E2m else None
    print(f"   Psi_A(C_1) = {Psi}  ~ {float(Psi):.10e};  E_A l^(2m) = {float(E2m):.10e};  rel = {float(ratio) if ratio is not None else None:.6e}"
          f"  ln Psi = {math.log(float(Psi)) if Psi > 0 else float('-inf'):.4f}  ({time.time()-t1:.0f}s)", flush=True)
print(f"done ({time.time()-t0:.0f}s)")
