#!/usr/bin/env python3
"""Referee: own complex128 isotypic engine (EVIDENCE: floating point).
A from enum2 (.bin int8 rows) or .npy. Eigenspace of character zeta^c (f(gk) = zeta^c f(k), (gk)_t = k_{t-1}).
Columns: P_c(k^alpha) on translation-orbit reps of A, alpha over monomial orbit reps with |alpha| < m.
Weighted inner product (weight q per free orbit, 1 for k=0). Residual of F_i = khat_1^i khat_{-1}^{m-i}
by SVD projection with relative tolerances tol in {1e-8, 1e-10, 1e-12}. Prints rank per eigenspace and Psi_A(C_1)."""
import sys, os, time, itertools
import numpy as np
from math import comb

def load(q, L):
    p = f'data/A_q{q}_L{L}'
    if os.path.exists(p + '.npy'): return np.load(p + '.npy').astype(np.int8)
    return np.fromfile(p + '.bin', dtype=np.int8).reshape(-1, q)

def orbit_reps(A, q):
    # canonical rep = lexicographically largest rotation, encoded base 3
    w = 3 ** np.arange(q - 1, -1, -1, dtype=np.int64)
    codes = np.stack([(np.roll(A, s, axis=1).astype(np.int64) + 1) @ w for s in range(q)], 1)
    own = (A.astype(np.int64) + 1) @ w
    isrep = own == codes.max(1)
    reps = A[isrep]
    nz = np.any(reps != 0, axis=1)
    return reps[nz], bool((~nz).any())

def mono_orbit_reps(q, m):
    reps = set()
    def rec(x, left, cur):
        if x == q:
            t = tuple(cur)
            best = max(t[s:] + t[:s] for s in range(q))
            reps.add(best); return
        for e in range(0, min(2, left) + 1):
            cur.append(e); rec(x + 1, left - e, cur); cur.pop()
    # enumerate by support pattern to stay cheap: alpha with |alpha| < m
    out = set()
    def gen(deg):
        # choose multiset of positions: n2 doubles, n1 singles, n1 + 2 n2 = deg
        for n2 in range(deg // 2 + 1):
            n1 = deg - 2 * n2
            for S2 in itertools.combinations(range(q), n2):
                rest = [x for x in range(q) if x not in S2]
                for S1 in itertools.combinations(rest, n1):
                    a = [0] * q
                    for x in S2: a[x] = 2
                    for x in S1: a[x] = 1
                    if a[0] == 0 and deg > 0:
                        continue  # canonical choice below handles; require 0 in support to cut by ~q
                    t = tuple(a)
                    out.add(max(t[s:] + t[:s] for s in range(q)))
    for deg in range(m):
        gen(deg)
    return sorted(out)

def run(q, L, ms, tols=(1e-8, 1e-10, 1e-12)):
    t0 = time.time()
    A = load(q, L); U = A.shape[0]
    R, has0 = orbit_reps(A, q)
    assert len(R) * q + (1 if has0 else 0) == U, (len(R), U)
    print(f'q={q} L={L}: U={U} free orbits={len(R)} ({time.time()-t0:.0f}s)', flush=True)
    zeta = np.exp(2j * np.pi / q)
    Rf = R.astype(np.float64)
    shifts = [np.roll(Rf, s, axis=1) for s in range(q)]   # (g^s r)_t = r_{t-s}
    kh1 = Rf @ zeta ** np.arange(q); khm = np.conj(kh1)
    for m in ms:
        mons = mono_orbit_reps(q, m)
        cs = sorted(set((2 * i - m) % q for i in range(m + 1)))
        cols = {c: np.zeros((len(R), len(mons)), dtype=np.complex128) for c in cs}
        for s in range(q):
            Xs = shifts[s]
            V = np.ones((len(R), len(mons)))
            for j, al in enumerate(mons):
                v = np.ones(len(R))
                for x, e in enumerate(al):
                    if e: v = v * Xs[:, x] ** e
                V[:, j] = v
            for c in cs:
                cols[c] += np.exp(-2j*np.pi*((c*s) % q)/q) * V / q
        res = {}
        for c in cs:
            Mc = cols[c] * np.sqrt(q)
            nrm = np.linalg.norm(Mc, axis=0); keep = nrm > 1e-8 * np.sqrt(Mc.shape[0])
            Mc = Mc[:, keep] / nrm[keep]
            if c == 0:   # add the fixed point k = 0 (weight 1): monomial values at 0 are 1 for alpha = 0 only
                row0 = np.array([1.0 if sum(al) == 0 else 0.0 for al in mons])[keep] / nrm[keep]
                Mc = np.vstack([Mc, row0[None, :]])
            Uu, sv, _ = np.linalg.svd(Mc, full_matrices=False)
            res[c] = (Uu, sv, Mc.shape[1])
        out = {}
        for tol in tols:
            Rn = []
            ranks = {}
            for i in range(m + 1):
                c = (2 * i - m) % q
                F = (kh1 ** i) * (khm ** (m - i)) * np.sqrt(q)
                if c == 0:
                    F = np.concatenate([F, [0.0 if m > 0 else 1.0]])
                Uu, sv, nc = res[c]
                r = int((sv > tol * sv[0]).sum()); ranks[c] = r
                Q = Uu[:, :r]
                resid = F - Q @ (Q.conj().T @ F)
                Rn.append(np.vdot(resid, resid).real / U)
            Psi = sum(comb(m, i) ** 2 * Rn[i] for i in range(m + 1)) / 4 ** m
            out[tol] = (Psi, ranks, Rn)
            print(f'  m={m} tol={tol:g}: ranks {ranks}  ||R_i||^2/U = ' + ', '.join(f'{x:.6e}' for x in Rn) +
                  f'  Psi_A(C_1) = {Psi:.10e}  ({time.time()-t0:.0f}s)', flush=True)

if __name__ == '__main__':
    q, L = int(sys.argv[1]), int(sys.argv[2])
    run(q, L, [int(x) for x in sys.argv[3].split(',')])
