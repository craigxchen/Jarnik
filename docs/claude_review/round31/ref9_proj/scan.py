"""Scan all circles x^2+y^2=N, N<=NMAX; maximal windows of arc length <= C sqrt(R) with M>=MMIN points.
For each window (deduplicated under the 8 lattice symmetries): normalise, exponent rows, and
 (i) the exact larger-sieve chain for b in a set B (small directions + certificate b's),
 (ii) certificates c = e_i+e_j, 2e_i, and (for M<=MSMALL) all c with entries in {-1,..,2}, ||c||_1<=CMAX, sum 2,
 (iii) the angle claim and the derivative bound for X_b.
Referee's own code."""
import sys; sys.path.insert(0, '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref9_proj')
from cert import *
import numpy as np, json, time
from math import isqrt, exp
NMAX = int(sys.argv[1]); C = float(sys.argv[2]); MMIN = int(sys.argv[3]); out = sys.argv[4]
MSMALL = 9; CMAX = 4
L = isqrt(NMAX)
xs_ = np.arange(-L, L+1, dtype=np.int64)
X, Y = np.meshgrid(xs_, xs_)
X = X.ravel(); Y = Y.ravel(); Nn = X*X+Y*Y
m = (Nn <= NMAX) & (Nn > 0)
X, Y, Nn = X[m], Y[m], Nn[m]
A = np.arctan2(Y, X)
order = np.lexsort((A, Nn))
X, Y, Nn, A = X[order], Y[order], Nn[order], A[order]
starts = np.flatnonzero(np.r_[True, Nn[1:] != Nn[:-1]])
ends = np.r_[starts[1:], len(Nn)]
SYM = [lambda z: z, lambda z: (-z[1], z[0]), lambda z: (-z[0], -z[1]), lambda z: (z[1], -z[0])]
def canon(cl):
    best = None
    for s in range(4):
        for cj in (False, True):
            t = []
            for z in cl:
                w = z
                for _ in range(s): w = (-w[1], w[0])
                if cj: w = (w[0], -w[1])
                t.append(w)
            t = tuple(sorted(t))
            if best is None or t < best: best = t
    return best
SMALLB = [(x, y) for x in range(0, 8) for y in range(0, 8) if (x, y) != (0, 0) and x*x+y*y <= 50 and gcd(x, y) == 1]
seen = set(); res = []
t0 = time.time()
for s, e_ in zip(starts, ends):
    k = e_ - s
    if k < MMIN: continue
    N = int(Nn[s]); R = sqrt(N); lim = C/sqrt(R)
    ang = A[s:e_]; pts = list(zip(X[s:e_].tolist(), Y[s:e_].tolist()))
    ang2 = np.r_[ang, ang + 2*pi]
    prevj = -1
    for i in range(k):
        j = int(np.searchsorted(ang2, ang2[i] + lim, side='right')) - 1
        j = min(j, i + k - 1)
        if j - i + 1 >= MMIN and j > prevj:
            cl = [pts[t % k] for t in range(i, j+1)]
            key = (N, canon(cl))
            if key not in seen:
                seen.add(key); res.append((N, cl))
        prevj = max(prevj, j)
print("windows", len(res), "time", time.time()-t0, flush=True)

stats = dict(windows=len(res), sieve_fail=0, angle_fail=0, deriv_fail=0, ident_fail=0,
             min_ratio_proj=1e9, min_ratio_proj_arg=None, min_ratio_proj_sharp=1e9,
             min_ratio_cert=1e9, min_ratio_cert_arg=None, rows=[])
for N, cl in res:
    g, zs, Np = normalise(cl)
    M = len(zs); R1 = sqrt(Np)
    a, th, wid = angles(zs)
    primes, pis, e, rows = rows_of(zs, Np) if Np > 1 else ([], [], [], [[] for _ in zs])
    s = len(primes); Pp = max(primes) if primes else 1
    cs = []
    for i, j in combinations(range(M), 2):
        c = [0]*M; c[i] += 1; c[j] += 1; cs.append(('pair', c))
    for i in range(M):
        c = [0]*M; c[i] = 2; cs.append(('single', c))
    if M <= MSMALL:
        for c in product(range(-1, 3), repeat=M):
            if sum(c) == 2 and sum(map(abs, c)) <= CMAX and sum(map(abs, c)) >= 3:
                cs.append(('small', list(c)))
    minpairL1 = None; minpairT = None; minT = None; minCert = None; best_c = None
    Bset = set(SMALLB)
    for kind, c in cs:
        if Np == 1: break
        T, b, u, v = certificate(zs, Np, primes, pis, e, rows, c)
        l1 = sum(map(abs, c))
        # angle claim: arg b = (1/2) sum c_i arg z_i mod pi, within l1*wid/4 of theta
        ab = atan2(b[1], b[0])
        half = 0.5*sum(ci*ai for ci, ai in zip(c, a))
        d1 = abs((ab - half + pi/2) % pi - pi/2)
        d2 = abs((ab - th + pi/2) % pi - pi/2)
        if d1 > 1e-6 or d2 > l1*wid/4 + 1e-9: stats['angle_fail'] += 1
        xs = projections(zs, b); Xb = xs[-1]-xs[0]
        bound = sqrt(nrm(b))*R1*wid*(abs(sin(th-ab)) + wid/2)
        if Xb > bound + 1e-6: stats['deriv_fail'] += 1
        Bset.add(b)
        val_ = (1+l1)*sqrt(T)/M**1.5
        if minCert is None or val_ < minCert: minCert = val_; best_c = (kind, c, T)
        if kind == 'pair' and (minpairT is None or T < minpairT): minpairT = T
        if kind == 'pair':
            L1v = sum(map(abs, v)); minpairL1 = L1v if minpairL1 is None else min(minpairL1, L1v)
        if minT is None or T < minT: minT = T
    if minCert is not None and minCert < stats['min_ratio_cert']:
        stats['min_ratio_cert'] = minCert; stats['min_ratio_cert_arg'] = (N, M, best_c)
    # projection chain
    bs = sum(log(p)/p for p in primes if p <= M)
    best_proj = None
    for b in Bset:
        xs = projections(zs, b)
        n = len(xs)
        if 2*n < M: stats['sieve_fail'] += 1
        if n >= 2:
            lhs, mid, rhs = sieve_terms(xs, Np*nrm(b), M)
            if not (lhs <= mid + 1e-9 and mid <= rhs + 1e-9): stats['sieve_fail'] += 1
        Xb = xs[-1]-xs[0]
        r = (2+Xb)/M**1.5
        if best_proj is None or r < best_proj[0]: best_proj = (r, b, Xb, (2+Xb)/(M**2*exp(-bs)))
    if best_proj[0] < stats['min_ratio_proj']:
        stats['min_ratio_proj'] = best_proj[0]; stats['min_ratio_proj_arg'] = (N, M, best_proj[1], best_proj[2])
    stats['min_ratio_proj_sharp'] = min(stats['min_ratio_proj_sharp'], best_proj[3])
    stats['rows'].append(dict(N=N, Np=Np, M=M, s=s, Pp=Pp, minpairT=minpairT, minpairL1=minpairL1, minT=minT, minCert=minCert,
                              bestproj=best_proj[0], bestb=best_proj[1], Xb=best_proj[2], arcC=R*wid/sqrt(R)))
print("time", time.time()-t0)
summ = {k: v for k, v in stats.items() if k != 'rows'}
print(json.dumps(summ, default=str))
json.dump(stats, open(out, 'w'), default=str)
