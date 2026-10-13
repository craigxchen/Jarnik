"""Referee: exact infimum of tau*r over equal-weight Plotkin-feasible (r,tau), tau >= w0, Prop 4.1.
For each r: tau_r = least tau with sum_{x<y} ceil(r/2 + D_xy/tau) <= r floor(M/2) ceil(M/2) (closed set [tau_r,inf)).
LB = min_r r*max(w0, tau_r). Compare with the author's eqlb_out.json."""
import json, math, sys
import numpy as np
from sympy import nextprime
DATA = '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round7/p3data_checks/data'

def w0_of(X):
    p = X + 1
    while True:
        p = nextprime(p - 1) if p > 2 else 2
        if p > X and p % 4 == 1: return math.log(p)
        p += 1

def Tsum(d, r, tau):
    return int(np.sum(np.ceil(r/2 + d/tau - 1e-12)))

def Tinf(d, r):
    return int(np.sum(np.where(d > 0, r//2 + 1, math.ceil(r/2))))

def lb_exact(D, w0):
    M = D.shape[0]; iu = np.triu_indices(M, 1); d = D[iu]
    cap = lambda r: r*(M//2)*((M+1)//2)
    # r_P at tau = w0
    r = 1
    while Tsum(d, r, w0) > cap(r): r += 1
    r_w0 = r
    best = r_w0*w0; arg = (r_w0, w0)
    for r in range(1, r_w0):
        if Tinf(d, r) > cap(r): continue
        # breakpoints: tau where r/2 + d/tau hits an integer n, n > r/2 : tau = d/(n - r/2)
        # find least tau >= w0 at which Tsum <= cap via bisection on sorted candidate set
        cands = set()
        for dv in d[d > 0]:
            n0 = math.floor(r/2) + 1
            # tau = dv/(n - r/2) >= w0 -> n - r/2 <= dv/w0
            nmax = math.floor(r/2 + dv/w0)
            for n in range(n0, nmax+1):
                cands.add(dv/(n - r/2))
        cands = sorted(c for c in cands if c >= w0)
        if Tsum(d, r, w0) <= cap(r):
            tau_r = w0
        else:
            lo, hi = 0, len(cands)-1
            # Tsum nonincreasing in tau; find least cand with Tsum<=cap (largest cand -> passes since beyond all breakpoints = Tinf)
            if not cands:
                continue
            if Tsum(d, r, cands[-1]) > cap(r):
                continue
            while lo < hi:
                mid = (lo+hi)//2
                if Tsum(d, r, cands[mid]) <= cap(r): hi = mid
                else: lo = mid+1
            tau_r = cands[lo]
        if r*tau_r < best: best = r*tau_r; arg = (r, tau_r)
    return best, arg

if __name__ == '__main__':
    auth = json.load(open('/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round7/p3data_checks/r2/eqlb_out.json'))
    print(type(auth), (list(auth[0].keys()) if isinstance(auth, list) else list(auth.keys())[:3]))
    A = {}
    for R in (auth if isinstance(auth, list) else auth.values()):
        A[(R['M'], R['X'], R['t'])] = R
    res = {}
    for M in [8, 12, 20, 24, 32]:
        for Xf in (2, 4):
            X = Xf*M; w0 = w0_of(X)
            ratios = []
            for t in range(12):
                D = np.load(f'{DATA}/D_construct_M{M}_t{t}_u1_X{Xf}.npy')
                lam2 = 2*np.linalg.eigvalsh(D)[-1]
                lb, arg = lb_exact(D, w0)
                a = A.get((M, X, t))
                alb = a.get('LB_eq', a.get('lb')) if a else None
                ratios.append(lb/lam2)
                if a is None or alb is None or not (alb <= lb + 1e-6 and alb >= lb*(1-0.011)):
                    print('CHECK', M, X, t, lb, arg, alb, flush=True)
            res[(M, X)] = (np.median(ratios), min(ratios), max(ratios))
            print(M, X, 'w0=%.4f' % w0, 'LB/2lam median %.3f [%.3f, %.3f]' % res[(M, X)], flush=True)
