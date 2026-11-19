"""Referee reproduction of Part V.1 (own arithmetic).  usage: python3 check_V1.py r [gcd] [factor]"""
import sys, itertools
sys.path.insert(0, __import__('os').path.dirname(__file__))
from f4arith import *

def par(x): return bin(x).count('1') & 1

def run(r, do_gcd, do_factor):
    m = 2 ** r; d = r + 1; Q = 2 ** d
    V = list(range(Q))
    forms = [w for w in range(1, Q) if w & 1]
    assert len(forms) == m
    H = {w: frozenset(v for v in V if par(v & w) == 0) for w in forms}
    for w in forms:
        assert 1 not in H[w] and len(H[w]) == Q // 2
        assert set(V) - H[w] == {v ^ 1 for v in H[w]}       # 1+H = V\H
    # S_J by enumeration
    def S(J): return sum(deg2(1 ^ v) for v in J)
    Smax = 0; pairs = list(itertools.combinations(forms, 2))
    Sdist = {}
    for a, b in pairs:
        J = H[a] & H[b]; assert len(J) == Q // 4
        s = S(J); Sdist[s] = Sdist.get(s, 0) + 1; Smax = max(Smax, s)
    assert Smax == (d - 2) * 2 ** (d - 2), (Smax, (d-2)*2**(d-2))
    L = max(d, 2 * Smax)
    assert L == (d - 2) * 2 ** (d - 1) == m * (r - 1)
    D = (m - 1) * L
    T = 1 << L
    def Yv(v): return (T ^ v, 1)
    def prod4(vs):
        z = (1, 0)
        for v in vs: z = mul4(z, Yv(v))
        return z
    pts = {w: prod4(sorted(H[w] - {0})) for w in forms}
    # N independently in F_2[t]: prod over even v != 0 of x^2+x+1 at x = t^L + v
    N = 1
    for v in V:
        if v & 1 == 0 and v != 0:
            x = T ^ v
            N = clmul(N, clmul(x, x) ^ x ^ 1)
    assert deg2(N) == (Q - 2) * L == 2 * D
    for w, z in pts.items():
        assert deg4(z) == D
        assert (z[0] >> D) & 1 == 1 and (z[1] >> D) & 1 == 0     # monic
        assert norm4(z) == N
        assert mul4(z, conj4(z)) == (N, 0)
        # conjugate = product over complementary coset minus {1}
        assert conj4(z) == prod4(sorted((set(V) - H[w]) - {1}))
    assert len(set(pts.values())) == m
    conjset = {conj4(z) for z in pts.values()}
    nconj_in = len(conjset & set(pts.values()))
    maxdiff = -1; mind = 10**9
    for a, b in pairs:
        J = H[a] & H[b]
        diff = add4(pts[a], pts[b])
        dd = deg4(diff)
        PJ1 = 1
        for v in J: PJ1 = clmul(PJ1, 1 ^ v)
        rhs = prod4(sorted(J - {0}))
        rhs = (clmul(rhs[0], PJ1), clmul(rhs[1], PJ1))
        assert diff == rhs                         # identity (V.1)
        assert dd == (Q // 4 - 1) * L + S(J)
        maxdiff = max(maxdiff, dd); mind = min(mind, dd)
    assert 2 * maxdiff == D
    # conjugates vs points: distance
    cd = min(deg4(add4(conj4(pts[a]), pts[b])) for a in forms for b in forms)
    print(f"m={m} d={d} L={L} D={D} parent={2*D} pairs={len(pairs)} maxdiff={maxdiff} (D/2={D//2}) "
          f"mindiff={mind} S-distribution={sorted(Sdist.items())} conj-in-cluster={nconj_in} "
          f"min deg(conj z - z')={cd}")
    if do_gcd:
        ws = list(forms)
        g = lgcd(to_list(pts[ws[0]]), to_list(pts[ws[1]]))
        for w in ws[2:]:
            g = lgcd(g, lmod(to_list(pts[w]), g)) if len(g) > 1 else g
        dg = len(g) - 1
        print(f"   gcd of all points: degree {dg}")
        if dg > 0:
            # factor g's norm over F_2 to identify it
            gl = from_list(g)
            print("   factorisation of Norm(g) over F_2:",
                  sorted((deg2(p), e) for p, e in factor2(norm4(gl)).items()))
        Dp = D - dg; mp = maxdiff - dg
        print(f"   after division: D'={Dp} maxdiff'={mp} c'={mp - Dp/2}")
        # also gcd of each pair
        pg = []
        for a, b in pairs[:min(len(pairs), 40)]:
            pg.append(len(lgcd(to_list(pts[a]), to_list(pts[b]))) - 1)
        print(f"   pairwise gcd degrees (first {len(pg)} pairs): min {min(pg)} max {max(pg)}")
    if do_factor:
        fac = factor2(N)
        chk = 1
        for p, e in fac.items():
            for _ in range(e): chk = clmul(chk, p)
        assert chk == N
        summ = {}
        bad = 0
        for p, e in fac.items():
            k = (deg2(p), e); summ[k] = summ.get(k, 0) + 1
            if deg2(p) % 2 == 1 and e % 2 == 1: bad += 1
        print(f"   N factorisation (degree, multiplicity): count = {sorted(summ.items())}")
        print(f"   odd-degree primes to odd power: {bad};  N square? {all(e%2==0 for e in fac.values())}")
        sp = [p for p in fac if deg2(p) % 2 == 0]
        print(f"   #split primes {len(sp)}, total split degree with mult {sum(deg2(p)*fac[p] for p in sp)}, "
              f"inert degree with mult {sum(deg2(p)*fac[p] for p in fac if deg2(p)%2)}")
        return fac, pts, N, L
    return None

if __name__ == '__main__':
    r = int(sys.argv[1])
    run(r, 'gcd' in sys.argv, 'factor' in sys.argv)
