"""FLOATING-POINT (complex128) isotypic computation for larger sdc_q boxes (EVIDENCE; validated against exact_iso.py
on sdc7 / sdc11).  For each translation character c in 0..(q-1)/2 (c and -c are conjugate) build, degree block by
degree block, a weighted-orthonormal basis of Poly_{<=e,c} from the columns P_c(k^alpha) (monomial orbit reps, exponents
<= 2), with block CGS2 + SVD rank truncation (rel. tol TOL); record dims and, for every m <= mmax (m < q):
  * ||R_i||^2 for the F_i = khat_1^i khat_-1^{m-i} landing in c (and its conjugate partner m-i in -c),
  * the residual vectors R_i (for the Chebyshev bracket on the circle),
  * the residuals of the c-components of l_u^m for NRAND random u (|u|_inf = 1).
Output: npz with everything + a text log.
usage: python3 float_iso.py q L mmax [nrand] [tol]
"""
import sys, time, itertools, math
import numpy as np

q, L, mmax = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
NRAND = int(sys.argv[4]) if len(sys.argv) > 4 else 6
TOL = float(sys.argv[5]) if len(sys.argv) > 5 else 1e-9
COMPS = [int(x) for x in sys.argv[6].split(',')] if len(sys.argv) > 6 else None   # optional subset of components (then only m = mmax assembled)
mmax = min(mmax, q - 1)
t0 = time.time()
reps = np.fromfile(f"data/A_q{q}_L{L}.bin", dtype=np.int8).reshape(-1, q)
fx = [r for r in reps if (r == r[0]).all()]
free = np.array([r for r in reps if not (r == r[0]).all()], dtype=np.int8)
assert len(fx) == 1 and not fx[0].any()
R = len(free); U = q * R + 1
smaxA = int((free != 0).sum(axis=1).max())
zeta = np.exp(2j * np.pi / q)
orb = np.stack([np.roll(free, s, axis=1) for s in range(q)])        # q x R x q  (g^s rep)
print(f"[float_iso] sdc{q} L={L}: |A|={U}, free orbits {R}, max support {smaxA}, mmax={mmax}, tol={TOL}", flush=True)


def mono_reps(deg):
    out = []; seen = set()
    for sz in range(0, min(deg, smaxA) + 1):
        if 2 * sz < deg: continue
        for S in itertools.combinations(range(q), sz):
            if sz and S[0] != 0: continue
            for ex in itertools.product((1, 2), repeat=sz):
                if sum(ex) != deg: continue
                a = [0] * q
                for x, e in zip(S, ex): a[x] = e
                a = tuple(a)
                if a in seen: continue
                for s in range(q): seen.add(a[-s:] + a[:-s] if s else a)
                out.append(a)
    return out


MONO = [mono_reps(e) for e in range(mmax)]
print("monomial orbit reps by degree:", [len(x) for x in MONO], f"({time.time()-t0:.0f}s)", flush=True)

kh1 = free.astype(np.float64) @ zeta ** np.arange(q)                 # khat_1 at reps
khm = np.conj(kh1)
rng = np.random.default_rng(12345)
UR = rng.uniform(-1, 1, (NRAND, q)); UR /= np.abs(UR).max(axis=1, keepdims=True)
# components of l_u^m at reps: g_c(rep) = (1/q) sum_s zeta^{-cs} g(g^s rep)
LU = np.einsum('srx,nx->nsr', orb.astype(np.float64), UR)            # NRAND x q x R


# orb2[..., x] = k_x, orb2[..., q + x] = |k_x| (= k_x^2), orb2[..., 2q] = 1 : monomial value = product of <= deg gathers
orb2 = np.concatenate([orb, np.abs(orb), np.ones(orb.shape[:2] + (1,), dtype=np.int8)], axis=2)


def block_cols(c, e, chunk=None):
    """yield column chunks P_c(k^alpha) at reps for monomial reps of degree e (complex, R x n)."""
    w = zeta ** (-c * np.arange(q))
    al = MONO[e]
    if chunk is None: chunk = max(64, int(4e8 / (q * R * max(1, e))))
    for st in range(0, len(al), chunk):
        A_ = np.array(al[st:st + chunk], dtype=np.int64).reshape(-1, q)
        n = len(A_); smx = max(1, int((A_ > 0).sum(axis=1).max()))
        I = np.full((n, smx), 2 * q, dtype=np.int64)
        for t, a in enumerate(A_):
            nz = np.nonzero(a)[0]
            I[t, :len(nz)] = np.where(a[nz] == 1, nz, q + nz)
        Mv = orb2[:, :, I[:, 0]].copy()
        for s_ in range(1, smx):
            Mv *= orb2[:, :, I[:, s_]]
        Mv = Mv.reshape(q, -1).astype(np.float64)
        out = (w.real @ Mv + 1j * (w.imag @ Mv)).reshape(R, n)
        yield out


results = {}
for c in (COMPS if COMPS is not None else range(0, (q - 1) // 2 + 1)):
    tc = time.time()
    extra = (c == 0)
    nrow = R + (1 if extra else 0)
    sw = np.full(nrow, math.sqrt(q));
    if extra: sw[-1] = 1.0
    Q = np.zeros((nrow, 0), dtype=complex)
    dims = []
    for m in range(1, mmax + 1):
        e = m - 1
        # add degree-e block
        for B in block_cols(c, e):
            if extra:
                B = np.vstack([B, np.full((1, B.shape[1]), q * 1.0 if e == 0 else 0.0)])
            B = B * sw[:, None]
            nb = np.linalg.norm(B, axis=0); keep = nb > 1e-8 * math.sqrt(nrow * q)   # exact zeros (e.g. P_c(1), c != 0) come out ~1e-15
            B = B[:, keep] / nb[keep]
            if B.shape[1] == 0: continue
            for _ in range(2):
                B = B - Q @ (Q.conj().T @ B)
            assert np.isfinite(B).all(), 'non-finite block'
            try:
                Uu, S, _ = np.linalg.svd(B, full_matrices=False)
            except np.linalg.LinAlgError:
                import scipy.linalg as _sla
                Uu, S, _ = _sla.svd(B, full_matrices=False, lapack_driver='gesvd')
            k = int((S > TOL).sum())          # absolute tolerance relative to unit-normalised columns
            if k:
                newv = Uu[:, :k]
                for _ in range(2):
                    newv = newv - Q @ (Q.conj().T @ newv)
                    newv, _r = np.linalg.qr(newv)
                Q = np.hstack([Q, newv])
            if Q.shape[1] >= nrow: break
        dims.append(Q.shape[1])
        # circle: i with 2i - m = c mod q
        i = ((m + c) * pow(2, -1, q)) % q
        ent = {}
        if i <= m:
            F = kh1 ** i * khm ** (m - i)
            if extra: F = np.append(F, 0.0)
            Fw = F * sw
            r = Fw - Q @ (Q.conj().T @ Fw)
            r = r - Q @ (Q.conj().T @ r)
            ent['i'] = i; ent['R2'] = float(np.vdot(r, r).real) / U; ent['F2'] = float(np.vdot(Fw, Fw).real) / U
            ent['Rvec'] = (r / sw)[:R]
        # random u: component c of l_u^m
        rr = []
        for n in range(NRAND):
            g = LU[n] ** m                                                # q x R
            gc = (zeta ** (-c * np.arange(q))) @ g / q                    # R
            if extra: gc = np.append(gc, 0.0 if m % 2 else 0.0)          # fixed point k=0: l_u^m(0) = 0
            gw = gc * sw
            r = gw - Q @ (Q.conj().T @ gw); r = r - Q @ (Q.conj().T @ r)
            rr.append((float(np.vdot(r, r).real) / U, float(np.vdot(gw, gw).real) / U))
        ent['rand'] = rr
        ent['dim'] = Q.shape[1]; ent['rows'] = nrow
        results[(c, m)] = ent
        if 'i' in ent: print(f'   c={c} m={m}: i={ent["i"]} ||R_i||^2={ent["R2"]:.10e} ||F_i||^2={ent["F2"]:.10e} dim={Q.shape[1]}', flush=True)
    print(f"c={c}: dims Poly_<m,c (m=1..{mmax}) = {dims} / rows {nrow}  ({time.time()-tc:.0f}s)", flush=True)
    del Q

# assemble
if COMPS is not None:
    print('component subset run: per-component values printed above'); sys.exit(0)
out = {}
for m in range(1, mmax + 1):
    R2 = {}; F2 = {}
    for c in range(0, (q - 1) // 2 + 1):
        ent = results[(c, m)]
        if 'i' in ent:
            i = ent['i']; R2[i] = ent['R2']; F2[i] = ent['F2']; R2[m - i] = ent['R2']; F2[m - i] = ent['F2']
    Psi = sum(math.comb(m, i) ** 2 * R2[i] for i in range(m + 1)) / 4 ** m
    E2 = sum(math.comb(m, i) ** 2 * F2[i] for i in range(m + 1)) / 4 ** m
    # random u: sum over components (c and -c equal for c != 0)
    pr = []
    for n in range(NRAND):
        tot = 0.0; totg = 0.0
        for c in range(0, (q - 1) // 2 + 1):
            a, b = results[(c, m)]['rand'][n]
            mult = 1 if c == 0 else 2
            tot += mult * a; totg += mult * b
        pr.append((tot, totg))
    dimP = sum((1 if c == 0 else 2) * results[(c, m)]['dim'] for c in range(0, (q - 1) // 2 + 1))
    out[m] = (Psi, E2, pr, dimP)
    print(f"m={m}: dim Poly_<m|_A = {dimP}, dim V = {U - dimP}; Psi_A(C_1) = {Psi:.10e}, E l^2m = {E2:.6e}, rel = {Psi/E2:.6e};"
          f" random u: Psi in [{min(p[0] for p in pr):.4e}, {max(p[0] for p in pr):.4e}], rel in "
          f"[{min(p[0]/p[1] for p in pr):.4e}, {max(p[0]/p[1] for p in pr):.4e}]", flush=True)
np.savez_compressed(f"out/float_q{q}_L{L}.npz", results=np.array([results], dtype=object), UR=UR, q=q, L=L)
print(f"done ({time.time()-t0:.0f}s)")
