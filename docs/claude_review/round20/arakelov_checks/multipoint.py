"""Multi-point contact data at finite places (evidence only).

For each arc cluster and each odd rational prime q not dividing N with q > M:
the q-adic reduction tree of the cluster is determined by pair valuations (ultrametric). We record
  T-contacts with |T| >= 3: subsets of >= 3 points in one residue class mod q (unit-1 collisions, i.e.
  z_i = z_j mod q in Z[i]; for split q this is the same as mod pi and mod pibar, two_thirds §3).
Checks (exact): ultrametricity of v_q(z_i - z_j) (min of two >= third) at every q.
Compares the observed number of (q, triple) incidences with the random-class expectation
  sum_q binom(M,3) / m_q^2 , m_q = q - 1 (split), q + 1 (inert)  over q in (M, Qmax].
Also the share of sum gamma_ij carried by q > M.
"""
import math
from itertools import combinations
from sympy import factorint, primerange
from gauss import sub, norm
from clusters import standard_clusters

def vq(q, z):
    # valuation of Gaussian integer z at the rational prime q (unit-1 collision depth): min over places above q
    if q % 4 == 3:
        a = 0
        x, y = z
        while x % q == 0 and y % q == 0:
            x //= q; y //= q; a += 1
        return a
    # split: collisions at pi and pibar coincide for points of equal norm prime to q; use v_q of both coords
    a = 0
    x, y = z
    while x % q == 0 and y % q == 0:
        x //= q; y //= q; a += 1
    return a

def main():
    out = open('multipoint_out.txt', 'w')
    out.write(('%-22s %3s %8s %8s %10s %10s %8s' % ('cluster', 'M', 'sumgam', 'gam_q>M', 'triples_obs', 'triples_exp', 'ultra')) + ' pairs_obs pairs_exp\n')
    for cl in standard_clusters():
        if cl['C'] is None:
            continue
        pts, N, M = cl['pts'], cl['N'], cl['M']
        fN = factorint(N)
        pairs = list(combinations(range(M), 2))
        qs = set()
        sumgam = 0.0
        for i, j in pairs:
            n = norm(sub(pts[i], pts[j]))
            f = factorint(n)
            qs |= set(f)
        # sum gamma via norms of reduced chords: use log Norm(z_i - z_j) - (W - d_ij)
        from gauss import gcd_g, exact_div
        gam_large = 0.0
        for i, j in pairs:
            c = exact_div(sub(pts[i], pts[j]), gcd_g(pts[i], pts[j]))
            fc = factorint(norm(c))
            sumgam += math.log(norm(c))
            gam_large += sum(e * math.log(q) for q, e in fc.items() if q > M and q not in fN)
        ultra = True
        triples = 0
        pairs_obs = 0
        for q in sorted(qs):
            if q == 2 or q in fN or q <= M:
                continue
            V = {}
            for i, j in pairs:
                V[(i, j)] = vq(q, sub(pts[i], pts[j]))
                pairs_obs += (V[(i, j)] >= 1)
            for i, j, k in combinations(range(M), 3):
                a, b, c = V[(i, j)], V[(j, k)], V[(i, k)]
                s = sorted([a, b, c])
                ultra &= (s[0] == s[1])          # ultrametric: the two smallest agree
                if min(a, b, c) >= 1:
                    triples += 1
        # random-class expectation, restricted (correctly) to q <= min chord of the triple:
        # q | (z_i - z_j) with q not | N means q | c_ij, forcing |c_ij| >= q: cap at the REDUCED chord.
        import bisect
        plist = [q for q in primerange(M + 1, 4 * 10**6) if q not in fN]
        pref = [0.0]
        for q in plist:
            m = q - 1 if q % 4 == 1 else q + 1
            pref.append(pref[-1] + 1.0 / m ** 2)
        exp = 0.0
        for i, j, k in combinations(range(M), 3):
            ch = min(math.sqrt(norm(exact_div(sub(pts[a], pts[b]), gcd_g(pts[a], pts[b])))) for a, b in ((i, j), (j, k), (i, k)))
            exp += pref[bisect.bisect_right(plist, ch)]
        pexp = 0.0
        for i, j in pairs:
            ch = math.sqrt(norm(exact_div(sub(pts[i], pts[j]), gcd_g(pts[i], pts[j]))))
            # first-moment expectation sum_q 1/m_q over M < q <= chord
            k = bisect.bisect_right(plist, ch)
            pexp += sum((1.0 / ((q - 1) if q % 4 == 1 else (q + 1))) for q in plist[:k])
        out.write('%-22s %3d %8.2f %8.2f %10d %10.3f %8s %6d %8.2f\n' % (cl['name'], M, sumgam, gam_large, triples, exp, ultra, pairs_obs, pexp))
        out.flush()
    out.close()
    print(open('multipoint_out.txt').read())

if __name__ == '__main__':
    main()
