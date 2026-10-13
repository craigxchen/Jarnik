"""Exact lower bound for EQUAL-WEIGHT profiles (all columns weight tau >= w0) of the pair system mu >= D.
For weight tau and r columns: d_xy >= T_xy = ceil(r/2 + D_xy/tau) and sum_{x<y} d_xy <= r floor(M/2) ceil(M/2).
r_P(tau) = least r passing this test (nonincreasing in tau).  On [a, b]: tau r_P(tau) >= a r_P(b).
LB_eq = min over a geometric grid w0 (1+h)^i, i < I, of a_i r_P(a_(i+1)), together with
tau_max r_P(inf) for tau >= tau_max.  D/tau thresholds use exact comparison D = 2 log m via floats with
margin 1e-9 (we round thresholds DOWN by subtracting 1e-9 before ceil, which can only lower the bound).
Output: per instance LB_eq / 2 lam_max, and the grid minimiser."""
import sys, math, json, glob, re
import numpy as np
sys.path.insert(0, '..')
from qcore import primes_1mod4_above

def rP(D, tau, iu, cap, r0=1):
    Dv = D[iu]
    for r in range(r0, 100000):
        T = np.ceil(r / 2.0 + Dv / tau - 1e-9)
        if T.sum() <= r * cap:
            return r

def rP_inf(D, iu, cap):
    Dv = D[iu]
    for r in range(1, 100000):
        T = np.where(Dv > 0, math.floor(r / 2) + 1, math.ceil(r / 2))
        if T.sum() <= r * cap:
            return r

def lb_eq(D, w0, h=0.01):
    M = D.shape[0]; iu = np.triu_indices(M, 1); cap = (M // 2) * ((M + 1) // 2)
    lam = float(np.linalg.eigvalsh(D)[-1])
    rinf = rP_inf(D, iu, cap)
    taus = [w0]
    while taus[-1] < 4 * D.max() + w0:
        taus.append(taus[-1] * (1 + h))
    rs = [rP(D, t, iu, cap) for t in taus]
    best = min(taus[i] * rs[i + 1] for i in range(len(taus) - 1))
    best = min(best, taus[-1] * rinf)
    i0 = int(np.argmin([taus[i] * rs[i] for i in range(len(taus))]))
    return best, lam, taus[i0] / w0, rs[i0], min(taus[i] * rs[i] for i in range(len(taus)))

if __name__ == '__main__':
    out = []
    for M in (8, 12, 20, 24, 32):
        for Xf in (2, 4):
            X = M * Xf; w0 = math.log(primes_1mod4_above(X, 1)[0])
            for t in range(12):
                D = np.load('../data/D_construct_M%d_t%d_u1_X%d.npy' % (M, t, Xf))
                lb, lam, f, r, grid = lb_eq(D, w0)
                out.append(dict(M=M, X=X, t=t, w0=w0, lam2=2 * lam, LB_eq=lb, grid_min=grid, f=f, r=r))
                print('M=%2d X=%3d t=%2d 2lam %.2f LB_eq %.2f ratio %.3f (grid min %.2f at tau=%.2f w0, r=%d)' % (
                    M, X, t, 2 * lam, lb, lb / (2 * lam), grid, f, r), flush=True)
    json.dump(out, open('eqlb_out.json', 'w'), indent=0)


def summary(path='eqlb_out.json'):
    """medians of LB_eq / 2 lam_max per (M, X) over the saved instances (printed by eqlb_summary)."""
    o = json.load(open(path))
    lines = []
    for M in sorted({d['M'] for d in o}):
        for X in sorted({d['X'] for d in o if d['M'] == M}):
            r = [d['LB_eq'] / d['lam2'] for d in o if d['M'] == M and d['X'] == X]
            f = [d['f'] for d in o if d['M'] == M and d['X'] == X]
            lines.append('M=%d X=%d n=%d LB_eq/2lam median %.3f [%.3f, %.3f]; tau*/w0 median %.2f' % (
                M, X, len(r), np.median(r), min(r), max(r), np.median(f)))
    return '\n'.join(lines)
