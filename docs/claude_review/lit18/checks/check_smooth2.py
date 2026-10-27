"""(C') smooth-number tuples vs the round-8 system (5.1): counts + sampled exact verification.
Run: python3 -I check_smooth2.py   (imports helpers from check_tools.py by exec of its header)"""
import math, random, itertools, sys
from fractions import Fraction
src = open(__file__.replace('check_smooth2.py', 'check_tools.py')).read().split("if __name__")[0]
exec(src)
random.seed(1)

def run(P, R, m, nsample=3000, ncells_proj=300):
    X = R * R
    W, s = s_integers(P, X)
    side = max(1, math.isqrt(R) // 2)
    cells = {}
    for w in W:
        if 4 * gnorm(w) < X:
            continue
        cells.setdefault((w[0] // side, w[1] // side), []).append(w)
    U = sum(len(v) for v in cells.values())
    ncell_tot = int(math.pi * (R * R - R * R / 4) / (side * side))
    Tm = sum(math.perm(len(v), m) for v in cells.values())
    Tm1 = sum(math.perm(len(v), m - 1) for v in cells.values())
    # sample tuples uniformly from cells weighted by perm count
    keys = [k for k, v in cells.items() if len(v) >= m]
    weights = [math.perm(len(cells[k]), m) for k in keys]
    dmin, dmax, bigg = float('inf'), 0.0, 0
    for _ in range(nsample):
        k = random.choices(keys, weights)[0]
        tup = random.sample(cells[k], m)
        g = tup[0]
        for w in tup[1:]:
            g = ggcd(g, w)
        if gnorm(g) > R:
            bigg += 1
        x = [gexact_div(w, g) for w in tup]
        H = max(gnorm(xi) for xi in x)
        fin = -sum(math.log(gnorm(xi)) for xi in x)
        arch = math.log(gnorm(x[0])) + sum(math.log(gnorm((xi[0]-x[0][0], xi[1]-x[0][1]))) for xi in x[1:]) - m*math.log(H)
        d = -(fin + arch) / math.log(H) - m
        dmin = min(dmin, d); dmax = max(dmax, d)
    # projective multiplicity in a sample of cells: ordered tuples vs distinct projective points
    tot = 0; proj = set()
    for k in random.sample(keys, min(ncells_proj, len(keys))):
        pts = cells[k]
        if math.perm(len(pts), m) > 20000:
            continue
        for tup in itertools.permutations(pts, m):
            tot += 1
            x0 = tup[0]
            proj.add(tuple((Fraction(gmul(xi, gconj(x0))[0], gnorm(x0)), Fraction(gmul(xi, gconj(x0))[1], gnorm(x0))) for xi in tup[1:]))
    print(f"P={P} s={s} R={R} m={m}: U={U} occupied cells={len(cells)} (of ~{ncell_tot}), mean mu={U/ncell_tot:.2f}; "
          f"ordered m-tuples={Tm}, ordered (m-1)-tuples={Tm1}, ratio Tm/Tm1={Tm/max(Tm1,1):.2f}; "
          f"sampled delta_eff in [{dmin:.3f},{dmax:.3f}] (window (1,{(m-1)/2}) for m=4), Norm(gcd)>R in {bigg}/{nsample}; "
          f"projective multiplicity in sampled cells = {tot/max(len(proj),1):.3f}")

for (P, R, m) in [(100, 10**4, 4), (200, 10**4, 4), (300, 3*10**4, 4), (300, 3*10**4, 5)]:
    run(P, R, m)
