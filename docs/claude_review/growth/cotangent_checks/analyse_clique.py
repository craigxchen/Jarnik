"""Exact anatomy of an integer cotangent clique: primitive points, cut profile,
residues, C_*, K5 = maxchord/N^(1/5), and the exponent log N/log(A/L)."""
import math, sys
from itertools import combinations
from cot_lib import (primitive_tuple, all_edge_lcm, factor, gnorm, ggcd, gsub,
                     split_prime_pi, gdivides, gexact, edge_norm, pair_cot)

def anatomy(xs, L):
    xs = sorted(xs)
    N = all_edge_lcm(xs, L)
    rows = primitive_tuple(xs, L)
    assert all(gnorm(z) == N for z in rows)
    A = xs[0]
    M = len(rows)
    Cs = 2 * N ** 0.25 * math.atan(L / A)
    chords = {}
    for i, j in combinations(range(M), 2):
        chords[(i, j)] = math.sqrt(gnorm(gsub(rows[i], rows[j])))
    K5 = max(chords.values()) / N ** 0.2
    print('L=%d X=%s' % (L, xs))
    print('N=%d = %s' % (N, factor(N)))
    print('C*=%.4f  maxchord/N^(1/5)=%.4f  log N/log(A/L)=%.4f' % (Cs, K5, math.log(N) / math.log(A / L)))
    # cut profile: for each split prime layer, the set of rows on the 'high' side
    prof = {}
    for p, e in factor(N).items():
        pi = split_prime_pi(p)
        al = []
        for z in rows:
            t = 0
            while gdivides(pi, z):
                z = gexact(z, pi); t += 1
            al.append(t)
        for tau in range(1, e + 1):
            side = frozenset(i for i in range(M) if al[i] >= tau)
            if len(side) > M / 2 or (len(side) == M / 2 and 0 not in side):
                side = frozenset(range(M)) - side
            prof[side] = prof.get(side, 0) + math.log(p)
    W = math.log(N)
    for side, w in sorted(prof.items(), key=lambda t: -t[1]):
        print('   cut %s|rest  weight/W=%.3f' % (sorted(side), w / W))
    res = []
    for i, j in combinations(range(M), 2):
        Q = xs[j - 1] if i == 0 else pair_cot(xs[i - 1], xs[j - 1], L)
        res.append(L // math.gcd(Q, L))
    print('residues s_ij:', res)
    return N

if __name__ == '__main__':
    L = int(sys.argv[1]); xs = [int(a) for a in sys.argv[2:]]
    anatomy(xs, L)
