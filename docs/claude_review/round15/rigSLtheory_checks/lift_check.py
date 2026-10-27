"""All-level one-modulus SL witnesses: for level-1 labels at q0 that are SL at phi_0 (classify.py DUMPSL), search a
lift ell2 mod m*q0 (ell2 = ell mod m) such that no v in V_K collides at level 2; then recompute the ALL-LEVEL demand
(genuine labels at all levels on F minus q0, ell / ell2 at q0; levels q0^a with a >= 3 cannot occur without level 2)
and verify exactly: PU, and (1.2) at phi = arg pi_j, i.e. m_(v,s) <= |D_s(Z_v)|; report whether the label lies outside
every eps-class and outside every (eps, t)-class (associate twists t, labels mod n are unaffected by t).
Usage: python3 lift_check.py FFILE q0 K dumpfile p1..pr"""
import sys, json, random, itertools, math
from sllib import two_sq, genuine, nu, Zv, Ds, VK, vq
from classify import load_F
Fpath, q0, K, dump = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]; ps = [int(x) for x in sys.argv[5:]]
F = load_F(Fpath); r = len(ps); pis = [two_sq(p) for p in ps]; V = VK(r, K)
m = q0-1 if q0 % 4 == 1 else q0+1; n = m//4; M2 = m*q0
_, L = genuine(ps, q0)
units = [u for u in range(1, n+1) if math.gcd(u, n) == 1] if n > 1 else [1]
epsc = set(tuple((u*e*x) % n for e, x in zip(eps, L)) for u in units for eps in itertools.product([1, -1], repeat=r))
pre = []
for v in V:
    z = Zv(pis, v); N = z[0]**2+z[1]**2; D = [Ds(z, s) for s in range(4)]
    g = []
    for s in range(4):
        gg = 1
        for q in F:
            if q != q0: gg *= q**vq(D[s], q)
        g.append(gg)
    pre.append((v, N, D, g))
ok = 0
for ell in json.load(open(dump)):
    assert [x % 2 for x in ell] == nu(ps, q0)
    random.seed(7); lift = None
    coll = [v for v in V if sum(a*b for a, b in zip(v, ell)) % n == 0]
    for _ in range(20000):
        cand = [x + m*random.randrange(q0) for x in ell]
        if all(sum(a*b for a, b in zip(v, cand)) % (n*q0) for v in coll): lift = cand; break
    if lift is None: print('  no lift found for', ell); continue
    pu = slab = 0
    for v, N, D, g in pre:
        d1 = sum(a*b for a, b in zip(v, ell)) % m
        for s in range(4):
            ms = g[s]*(q0 if d1 == s*n else 1)
            pu += ms*ms > N*(2 if s % 2 else 1); slab += ms > abs(D[s])
    out = tuple(x % n for x in ell) not in epsc
    good = pu == 0 and slab == 0 and out
    ok += good
    print(f"  ell={ell} lift={lift} outside every eps-class={out} PU viol={pu} slab viol at phi0={slab}")
print(f"ps={ps} K={K} q0={q0}: verified all-level SL witnesses outside every eps-class: {ok}")
