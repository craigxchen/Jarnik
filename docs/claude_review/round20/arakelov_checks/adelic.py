"""Theorem B checks: the (2/3+eps) finite inequality as an adelic M-point (Fekete) diameter inequality,
and the place-by-place slack of actual extremal clusters.

For each normalised cluster and each rational prime q:
  act(q)  = sum_{i<j} sum_{varpi | q} v_varpi(z_i - z_j) log Norm(varpi)          (actual local energy)
  min(q)  = exact local minimum of that energy over M-tuples of local points with the cluster's local data:
            q = 2      : lower bound binom(M,2) log 2 (only (1+i) | z_i - z_j is used, as in two_thirds)
            inert l    : sum_a E(M,(l+1) l^(a-1)) * 2 log l
            split q!|N : sum_a E(M,(q-1) q^(a-1)) * 2 log q
            p | N      : min over level vectors (n_0..n_e) of
                         [ e(M^2/4 - M/2) + sum_tau (S_tau - M/2)^2 + 2 sum_t sum_a E(n_t,(p-1)p^(a-1)) ] log p
  archimedean: act(inf) = - sum_{i<j} log|z_i - z_j|^2 ;
               min(inf) = - 2 binom(M,2) log delta_M(arc)   (Fekete points of the spanned chord, Lobatto nodes)
               crude    = - 2 binom(M,2) log(chord)          (the two_thirds bound uses the chord <= C sqrt R)
Checks: sum_v act(v) = 0 (product formula, to float precision); every act >= min (exact for finite q).
Output: slack decomposition. Excess at q > M (q not | N) is the 'large-prime part' of sum gamma_ij.
"""
import math
from itertools import combinations, product
from sympy import factorint, primerange
from gauss import sub, conj, norm, valuation, split_prime_factor
from clusters import standard_clusters
from dictionary import gaussian_primes_over
from numpy.polynomial import legendre as L
import numpy as np

def E(n, Q):
    if n <= 0:
        return 0
    b, r = divmod(n, Q)
    return Q * b * (b - 1) // 2 + r * b

def coll_min(n, base, p):
    """sum_a E(n, base * p^(a-1)) for a >= 1 until classes exceed n."""
    s, a = 0, 1
    while True:
        Q = base * p ** (a - 1)
        if Q >= n and E(n, Q) == 0:
            break
        s += E(n, Q)
        a += 1
    return s

def split_level_min(M, e, p):
    best = None
    def rec(prefix, remaining, slots):
        nonlocal best
        if slots == 1:
            ns = prefix + [remaining]
            S = 0; cost = 0.0
            for t in range(e):
                S += ns[t]
                cost += (S - M / 2) ** 2
            cost += 2 * sum(coll_min(x, p - 1, p) for x in ns)
            if best is None or cost < best:
                best = cost
            return
        for x in range(remaining + 1):
            rec(prefix + [x], remaining - x, slots - 1)
    rec([], M, e + 1)
    return e * (M * M / 4 - M / 2) + best

def lobatto_delta(M):
    """M-th diameter of [-1,1] (Fekete = Lobatto nodes), returned as delta_M; for [0,L] multiply by L/2."""
    if M == 2:
        return 2.0
    c = np.zeros(M); c[-1] = 1           # P_{M-1}
    dP = L.legder(c)
    x = np.concatenate(([-1.0], np.sort(L.legroots(dP)), [1.0]))
    s = sum(math.log(abs(x[i] - x[j])) for i in range(M) for j in range(i + 1, M))
    return math.exp(s / (M * (M - 1) / 2))

def analyse(cl):
    pts, N, M = cl['pts'], cl['N'], cl['M']
    fN = factorint(N)
    pairs = list(combinations(range(M), 2))
    primes = set([2]) | set(fN)
    for i, j in pairs:
        primes |= set(factorint(norm(sub(pts[i], pts[j]))))
    for q in primerange(3, M + 2):
        primes.add(q)
    rows = {}
    for q in sorted(primes):
        act = 0.0
        for varpi, Nv in gaussian_primes_over(q):
            for i, j in pairs:
                act += valuation(varpi, sub(pts[i], pts[j])) * math.log(Nv)
        if q == 2:
            mn = len(pairs) * math.log(2)
        elif q in fN:
            mn = split_level_min(M, fN[q], q) * math.log(q)
        elif q % 4 == 3:
            mn = coll_min(M, q + 1, q) * 2 * math.log(q)
        else:
            mn = coll_min(M, q - 1, q) * 2 * math.log(q)
        assert act >= mn - 1e-9, (cl['name'], q, act, mn)
        rows[q] = (act, mn)
    # archimedean
    act_inf = -sum(math.log(norm(sub(pts[i], pts[j]))) for i, j in pairs)
    angs = sorted(math.atan2(z[1], z[0]) for z in pts)
    # spanned chord (cluster sorted by angle within an arc < pi)
    a = [math.atan2(z[1], z[0]) for z in pts]
    a0 = min(a);
    span = max(((x - a0) % (2 * math.pi)) for x in a)
    if span > math.pi:            # wrap: recompute relative to max gap
        a_sorted = sorted(a); gaps = [(a_sorted[(k + 1) % M] - a_sorted[k]) % (2 * math.pi) for k in range(M)]
        span = 2 * math.pi - max(gaps)
    R = math.sqrt(N)
    chord = 2 * R * math.sin(span / 2)
    deltaM = (chord / 2) * lobatto_delta(M)
    min_inf = -2 * len(pairs) * math.log(deltaM)
    crude_inf = -2 * len(pairs) * math.log(chord)
    total_act = act_inf + sum(v[0] for v in rows.values())
    total_min_fek = min_inf + sum(v[1] for v in rows.values())
    ex_inf = act_inf - min_inf
    ex_small = sum(v[0] - v[1] for q, v in rows.items() if (q <= M or q == 2) and q not in fN)
    ex_N = sum(v[0] - v[1] for q, v in rows.items() if q in fN)
    ex_large = sum(v[0] - v[1] for q, v in rows.items() if q > M and q not in fN)
    return dict(name=cl['name'], M=M, W=math.log(N), C=cl['C'], prodformula=total_act,
                slack=-total_min_fek, ex_inf=ex_inf, fek_gain=min_inf - crude_inf,
                ex_small=ex_small, ex_N=ex_N, ex_large=ex_large,
                sum_gamma=sum(math.log(norm(sub(pts[i], pts[j]))) for i, j in pairs))

if __name__ == '__main__':
    with open('adelic_out.txt', 'w') as f:
        f.write('%-22s %3s %7s %7s %10s %8s %8s %8s %8s %8s %8s\n' % (
            'cluster', 'M', 'W', 'C', 'prodform', 'slack', 'ex_inf', 'fekgain', 'ex_q<=M', 'ex_p|N', 'ex_q>M'))
        for cl in standard_clusters():
            r = analyse(cl)
            f.write('%-22s %3d %7.2f %7s %10.2e %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n' % (
                r['name'], r['M'], r['W'], ('%.3f' % r['C']) if r['C'] else '-', r['prodformula'], r['slack'],
                r['ex_inf'], r['fek_gain'], r['ex_small'], r['ex_N'], r['ex_large']))
            f.flush()
        f.write('lobatto delta_M([-1,1]) for M=2..12, 50, 200: ' +
                ' '.join('%.5f' % lobatto_delta(m) for m in list(range(2, 13)) + [50, 200]) + '\n')
    print(open('adelic_out.txt').read())
