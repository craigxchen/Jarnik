"""Circle test (barrier.md §3.2-3.3) in toy sdc_q boxes: does the real circle
C_a = {u_theta(t) = cos(theta + 2 pi a t/q)} lie in TC(A, m)?  C_a in TC iff for every i = 0..m the function
F_i(k) = khat(a)^i * conj(khat(a))^(m-i),  khat(a) = sum_t zeta^(a t) k_t,
lies in Poly_{<m}|_A (x) C.  F_i is an eigenfunction of the translation with character zeta^(a(2i-m)); the test is done
in that eigenspace (rows: orbit representatives; for the trivial character the fixed points are added as rows;
columns: projected monomial orbit representatives of degree < m).  Arithmetic mod two primes p = 1 mod q (EVIDENCE for
char-0 ranks, exact mod p).  Also prints the deficiency delta_c = rows - rank for every character c (sum = dim V(A,m)).
Usage: python3 sdc_circles.py q L m [a]
"""
import sys, time
import numpy as np
import flint
from flib import sdc_patterns, box, reduced_monomials

q, L, m = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
a_ = int(sys.argv[4]) if len(sys.argv) > 4 else 1
t0 = time.time()
M, pats = sdc_patterns(q)
A = np.array(box(M, pats, L), dtype=np.int8)
Aset = {tuple(k): idx for idx, k in enumerate(A)}
seen = np.zeros(len(A), bool); reps = []; fixedpts = []
for idx in range(len(A)):
    if seen[idx]: continue
    orb = [Aset[tuple(np.roll(A[idx], s))] for s in range(q)]
    for j in orb: seen[j] = True
    (fixedpts if len(set(orb)) == 1 else reps).append(A[idx])
print(f"sdc{q} L={L} m={m}: |A|={len(A)} free orbits {len(reps)} fixed {len(fixedpts)} ({time.time()-t0:.0f}s)", flush=True)


def mono_reps(deg):
    out = []; seenm = set()
    for a in reduced_monomials(q, deg):
        if a in seenm: continue
        for s in range(q): seenm.add(tuple(np.roll(a, s)))
        out.append(np.array(a))
    return out


MONS = [mm for deg in range(m) for mm in mono_reps(deg)]
print(f"monomial orbit reps of degree < {m}: {len(MONS)}", flush=True)


def run(p):
    g0 = 2
    while pow(g0, (p - 1) // q, p) == 1: g0 += 1
    z = pow(g0, (p - 1) // q, p)
    defic = {}; member = {}
    need = sorted({(a_ * (2 * i - m)) % q for i in range(m + 1)})
    for c in need:
        rowsK = list(reps) + (list(fixedpts) if c == 0 else [])
        Rk = np.array(rowsK, dtype=np.int64)
        nr = len(Rk)
        chiinv = pow(pow(z, c, p), p - 2, p)
        shifts = [np.roll(Rk, s, axis=1) for s in range(q)]
        cols = []
        for al in MONS:
            acc = np.zeros(nr, dtype=np.int64)
            for s in range(q):
                v = np.prod(np.where(al > 0, shifts[s] ** al, 1), axis=1)    # monomial at g^s k
                acc = (acc + (v % p) * pow(chiinv, s, p)) % p
            cols.append(acc)
        # functions F_i with character c: a(2i - m) = c mod q
        fcols = {}
        zt = [pow(z, (a_ * t) % q, p) for t in range(q)]
        zbt = [pow(z, (-a_ * t) % q, p) for t in range(q)]
        for i in range(m + 1):
            if (a_ * (2 * i - m) - c) % q == 0:
                kh = [(sum(int(k[t]) * zt[t] for t in range(q))) % p for k in Rk]
                kb = [(sum(int(k[t]) * zbt[t] for t in range(q))) % p for k in Rk]
                fcols[i] = np.array([(pow(x, i, p) * pow(y, m - i, p)) % p for x, y in zip(kh, kb)], dtype=np.int64)
        Mt = flint.nmod_mat(nr, len(cols), [int(col[r]) for r in range(nr) for col in cols], p)
        rk = Mt.rank()
        defic[c] = nr - rk
        for i, f in fcols.items():
            Mf = flint.nmod_mat(nr, len(cols) + 1, [int(col[r]) for r in range(nr) for col in cols + [f]], p)
            member[i] = (Mf.rank() == rk)
    return defic, member


primes = [pp for pp in range(10 ** 9, 10 ** 9 + 400000) if (pp - 1) % q == 0 and flint.fmpz(pp).is_prime()][:2]
res = [run(p) for p in primes]
d1, m1 = res[0]; d2, m2 = res[1]
print("character deficiencies delta_c on the characters needed (c: delta) =", d1, "| second prime agrees:", d1 == d2)
print("membership of F_i (i: in span) =", m1, "| second prime agrees:", m1 == m2)
print(f"=> circle C_{a_} {'IS' if all(m1.values()) else 'is NOT'} contained in TC(A,{m}) (mod-p evidence)  "
      f"({time.time()-t0:.0f}s)")
