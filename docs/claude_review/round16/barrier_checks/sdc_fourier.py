"""Circle test of barrier.md §3.3 in small sdc_q boxes (toy models, NOT regime): is the complex Fourier direction
psi_1 = (zeta^t)_{t in F_q} in TC(A, i), i.e. is f_i(k) = (sum_t zeta^t k_t)^i in Poly_{<i}|_A ?
Equivalently: nu_A(psi_1) >= i iff f_i notin span.  Translation g: (g k)_t = k_{t-1}; f_i(g k) = zeta^i f_i(k), so the
test lives in the zeta^i-eigenspace of functions on A: rows = free g-orbits of A, columns = P_chi(monomial) over
monomial orbit representatives of degree < i (reduced monomials, exponents <= 2), plus the column f_i.
Arithmetic mod a prime p = 1 mod q (zeta in F_p).  Direction of validity: "f_i notin span mod p" implies
"notin span over Q(zeta)" provided rank_p(M) = rank_char0(M); we certify rank_p(M) = rank_p'(M) at a second prime as a
consistency check (EVIDENCE for the char-0 rank; exact otherwise).  Also reports rho_chi(i) := whether M is already
of full row rank (then the test is vacuous: TC = everything in that eigenspace).
Usage: python3 sdc_fourier.py q L imax
"""
import sys, time
import numpy as np
import flint
from flib import sdc_patterns, box, reduced_monomials

q, L, imax = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
t0 = time.time()
M, pats = sdc_patterns(q)
A = np.array(box(M, pats, L), dtype=np.int8)      # coordinates t = 0..q-1 (row t+1 = F_q element t)
print(f"sdc{q} L={L}: |A|={len(A)}  density={len(A)/3**q:.3e}  ({time.time()-t0:.0f}s)", flush=True)
# orbit representatives of A under cyclic shift
Aset = {tuple(k): idx for idx, k in enumerate(A)}
seen = np.zeros(len(A), bool); reps = []; fixed = 0
for idx in range(len(A)):
    if seen[idx]: continue
    k = A[idx]; orb = []
    for s in range(q):
        ks = tuple(np.roll(k, s)); j = Aset[ks]; orb.append(j)
    for j in orb: seen[j] = True
    if len(set(orb)) == 1: fixed += 1
    else:
        assert len(set(orb)) == q
        reps.append(k)
R = np.array(reps, dtype=np.int64)
print(f"free orbits {len(R)}, fixed points {fixed}", flush=True)
nzR = np.stack([((np.roll(R, s, axis=1) != 0) * (1 << np.arange(q))).sum(1) for s in range(q)])   # q x nrep
ngR = np.stack([((np.roll(R, s, axis=1) < 0) * (1 << np.arange(q))).sum(1) for s in range(q)])


def popparity(x):
    x = x.copy(); p = np.zeros_like(x)
    while np.any(x):
        p ^= (x & 1); x >>= 1
    return p


def mono_reps(deg):
    out = []; seenm = set()
    for a in reduced_monomials(q, deg):
        if a in seenm: continue
        orb = [tuple(np.roll(a, s)) for s in range(q)]
        for b in orb: seenm.add(b)
        out.append(a)
    return out


def test(i, p):
    F = flint.nmod
    # zeta: element of order q mod p
    g0 = 2
    while True:
        z = pow(g0, (p - 1) // q, p)
        if z != 1: break
        g0 += 1
    chi = pow(z, i, p); chiinv = pow(chi, p - 2, p)
    cols = []
    for deg in range(i):
        for a in mono_reps(deg):
            a = np.array(a)
            supp = int(((a > 0) * (1 << np.arange(q))).sum()); odd = int(((a % 2 == 1) * (1 << np.arange(q))).sum())
            col = np.zeros(len(R), dtype=object)
            acc = np.zeros(len(R), dtype=np.int64)
            for s in range(q):
                ok = (nzR[s] & supp) == supp
                sgn = 1 - 2 * popparity(ngR[s] & odd)
                v = np.where(ok, sgn, 0)               # monomial value at g^s k
                w = pow(chiinv, s, p)
                acc = (acc + (v % p) * w) % p
            cols.append(acc)
    # f_i on reps: (sum_t z^t k_t)^i
    zt = np.array([pow(z, t, p) for t in range(q)], dtype=object)
    fv = [(int(sum(int(R[r, t]) * zt[t] for t in range(q) if R[r, t])) % p) for r in range(len(R))]
    fv = np.array([pow(v, i, p) for v in fv], dtype=np.int64)
    Mt = flint.nmod_mat(len(R), len(cols), [int(c[r]) for r in range(len(R)) for c in cols], p) if cols else None
    rk = Mt.rank() if cols else 0
    Mf = flint.nmod_mat(len(R), len(cols) + 1, [int(c[r]) for r in range(len(R)) for c in cols + [fv]], p)
    rkf = Mf.rank()
    return len(cols), rk, rkf


for i in range(1, imax + 1):
    t1 = time.time()
    res = []
    for p in (1000000007 if (1000000007 - 1) % q == 0 else None, ):
        pass
    primes = [pp for pp in range(10 ** 9, 10 ** 9 + 200000) if (pp - 1) % q == 0 and flint.fmpz(pp).is_prime()][:2]
    out = [test(i, pp) for pp in primes]
    ncol, rk, rkf = out[0]
    verdict = "f_i NOT in span (nu_A(psi_1) >= i)" if rkf > rk else "f_i in span (psi_1 in TC(A,i))"
    full = " [M full row rank: vacuous]" if rk == len(R) else ""
    print(f"i={i}: monomial orbit reps {ncol}, rank_p {rk} / rows {len(R)}, rank with f_i {rkf}; second prime ranks "
          f"{out[1][1]}, {out[1][2]} -> {verdict}{full}  ({time.time()-t1:.0f}s)", flush=True)
