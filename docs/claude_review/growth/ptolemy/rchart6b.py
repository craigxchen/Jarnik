r"""Variant of rchart6.py: the leading coefficients of q_2..q_6 are FIXED (random, distinct,
nonzero) to exclude the degree-drop degeneration, and only 6 triple points are fixed.
Unknowns: a_1..a_6, r_2..r_6, lower coefficients of q_2..q_6 (5*9), 14 triple + 15 quad points
= 6+5+45+29 = 85 = number of equations.  NUMERICAL exploration only.
Usage: python3 rchart6b.py trials seed
"""
import sys, itertools
import numpy as np
import rchart6 as R

def build(rng):
    c = lambda n: (rng.normal(size=n) + 1j * rng.normal(size=n)) / np.sqrt(2)
    fixedT = R.choose_fixed(rng, nfix=6)
    fixedv = dict(zip(fixedT, c(len(fixedT)) * 2))
    freeB = [T for T in R.triples if T not in fixedT] + R.quads
    conds = []
    for T in R.triples + R.quads:
        Ts = sorted(T)
        for a, b in zip(Ts, Ts[1:]):
            conds.append((T, a, b))
    nq = R.DEG + 1
    leads = dict(zip(R.rows[1:], c(5)))
    S = dict(fixedv=fixedv, freeB=freeB, conds=conds, nq=nq, leads=leads)
    S['idx_a'] = {i: i - 1 for i in R.rows}
    S['idx_r'] = {i: 6 + (i - 2) for i in R.rows[1:]}
    S['idx_ql'] = {i: 11 + (i - 2) * (nq - 1) for i in R.rows[1:]}   # lower coefficients
    S['idx_p'] = {T: 11 + 5 * (nq - 1) + n for n, T in enumerate(freeB)}
    S['nz'] = 11 + 5 * (nq - 1) + len(freeB)
    assert S['nz'] == len(conds)
    S['z0'] = np.concatenate([c(6) * 2, c(5), c(5 * (nq - 1)) * 0.3, c(len(freeB)) * 2])
    return S

def full(S, z):
    """embed into rchart6 layout (with q including leads)"""
    nq = S['nq']
    zz = np.zeros(11 + 5 * nq + len(S['freeB']), dtype=complex)
    zz[:11] = z[:11]
    for i in R.rows[1:]:
        st = 11 + (i - 2) * nq
        zz[st] = S['leads'][i]
        zz[st + 1:st + nq] = z[S['idx_ql'][i]:S['idx_ql'][i] + nq - 1]
    zz[11 + 5 * nq:] = z[11 + 5 * (nq - 1):]
    T6 = dict(fixedv=S['fixedv'], freeB=S['freeB'], conds=S['conds'], nq=nq,
              idx_a=S['idx_a'], idx_r=S['idx_r'],
              idx_q={i: 11 + (i - 2) * nq for i in R.rows[1:]},
              idx_p={T: 11 + 5 * nq + n for n, T in enumerate(S['freeB'])}, nz=len(zz))
    return T6, zz

def F(S, z):
    T6, zz = full(S, z)
    return R.F(T6, zz)

def J(S, z):
    T6, zz = full(S, z)
    Jf = R.J(T6, zz)
    nq = S['nq']
    keep = list(range(11))
    for i in R.rows[1:]:
        st = 11 + (i - 2) * nq
        keep += list(range(st + 1, st + nq))
    keep += list(range(11 + 5 * nq, len(zz)))
    return Jf[:, keep]

def newton(S, z, iters=300):
    for it in range(iters):
        f = F(S, z); nf = np.linalg.norm(f)
        if nf < 1e-13:
            break
        step = np.linalg.lstsq(J(S, z), -f, rcond=None)[0]
        t = 1.0
        while t > 1e-9:
            zn = z + t * step
            if np.linalg.norm(F(S, zn)) < nf:
                break
            t /= 2
        z = zn
    return z, np.linalg.norm(F(S, z))

if __name__ == "__main__":
    trials = int(sys.argv[1]); seed = int(sys.argv[2])
    rng = np.random.default_rng(seed)
    good = 0
    for tr in range(trials):
        S = build(rng)
        z, nf = newton(S, S['z0'])
        if nf > 1e-9:
            print(f"trial {tr}: no convergence ({nf:.1e})", flush=True); continue
        T6, zz = full(S, z)
        np.save(f"r6bsol_seed{seed}_tr{tr}.npy", zz)
        out, why = R.analyse(T6, zz)
        if out is None:
            print(f"trial {tr}: converged; {why}", flush=True); continue
        sep, lsep, rmin, special, sep_big, rem = out
        ok = sep > 1e-5 and lsep > 1e-5 and rmin > 1e-5 and rem < 1e-8
        good += ok
        print(f"trial {tr}: converged; min sep {sep:.2e} (big {sep_big:.2e}), lead sep {lsep:.2e}, min|r| {rmin:.2e}, rem {rem:.1e} -> {'NONDEGENERATE' if ok else 'degenerate'}", flush=True)
    print("nondegenerate:", good, "/", trials)
