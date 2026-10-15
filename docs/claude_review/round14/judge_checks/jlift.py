"""Judge check J2b: complete three level-1 SL survivors at q0=107 (first4, K=4, from jonemod_out.txt) to all-level
systems: search a lift mod m*q0 with no level-2 collision, then recheck PU, ANG and the slab system at phi = arg pi_j
(m_(v,s) <= |D_s(Z_v)|) exactly with the all-level demand (genuine labels at all levels on F \\ {107})."""
import sys, os, itertools, random
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from jglib import two_sq, genuine, Zv, Ds, VK, vq, nu
from jang import box_meets_circle
ps = [5, 13, 17, 29]; K = 4; q0 = 107
F = eval([l for l in open('../referee_rigtheory_checks/forced_K4.txt') if l.startswith('F ')][0][2:])
pis = [two_sq(p) for p in ps]; V = VK(4, K); m = q0 + 1; n = m // 4; M2 = m * q0
_, L = genuine(ps, q0)
units = [u for u in range(1, n) if u % 3]
for ell in [(1, 4, 61, 92), (1, 88, 43, 50), (1, 88, 97, 104)]:
    assert [x % 2 for x in ell] == nu(ps, q0)
    outside = all(tuple((u*e*x) % n for e, x in zip(eps, L)) != tuple(x % n for x in ell)
                  for u in units for eps in itertools.product([1, -1], repeat=4))
    random.seed(1); lift = None
    for _ in range(5000):
        cand = [x + m*random.randrange(q0) for x in ell]
        if all(sum(a*b for a, b in zip(v, cand)) % (n*q0) for v in V): lift = cand; break
    if lift is None:  # collisions at level 1 force level-2 test only for colliding v
        for _ in range(20000):
            cand = [x + m*random.randrange(q0) for x in ell]
            coll = [v for v in V if sum(a*b for a, b in zip(v, ell)) % n == 0]
            if all(sum(a*b for a, b in zip(v, cand)) % (n*q0) for v in coll): lift = cand; break
    pu = ang = slab = 0
    for v in V:
        z = Zv(pis, v); N = z[0]**2 + z[1]**2; ms = []
        for s in range(4):
            d = Ds(z, s); g = 1
            for q in F:
                if q != q0: g *= q**vq(d, q)
            dd = sum(a*b for a, b in zip(v, lift)) % M2
            if dd % (n*q0) == 0: g *= q0**2  # cannot happen if lift is level-2 free
            elif sum(a*b for a, b in zip(v, ell)) % m == s*n: g *= q0
            ms.append(g)
            pu += g*g > N*(2 if s % 2 else 1); slab += g > abs(d)
        ang += not box_meets_circle(Fr(2*ms[2]**2, N)-1, 1-Fr(2*ms[0]**2, N), Fr(ms[3]**2, N)-1, 1-Fr(ms[1]**2, N))
    print('ell', ell, 'lift', lift, 'outside every eps-class', outside, '| PU viol', pu, 'ANG fail', ang,
          'slab viol at phi=arg pi', slab)
