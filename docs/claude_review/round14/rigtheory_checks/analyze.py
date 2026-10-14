"""Anatomy of non-genuine one-modulus survivors (round 14).  Usage: python3 analyze.py K "ps" q0 [nex]"""
import sys, math
from rlib14 import *
from forced import forced_set
from deform import run_q0
K = int(sys.argv[1]); ps = [int(x) for x in sys.argv[2].split(",")]; q0 = int(sys.argv[3])
nex = int(sys.argv[4]) if len(sys.argv) > 4 else 2
F = forced_set(ps, K, 60000)
d = run_q0(ps, K, F, q0, maxex=nex)
pis = gaussian_primes(ps); m, n, L = d['m'], d['n'], d['L']
print({k: d[k] for k in ['q0', 'n', 'total', 'inclass', 'strict', 'F2', 'F1']})
V = VK_mod_sign(len(ps), K)
others = [q for q in F if q != q0]
def coll(ell):
    out = []
    for v in V:
        t = sum(a * b for a, b in zip(v, ell)) % m
        if t % n == 0: out.append((v, t // n))
    return out
gen = coll(L)
print("genuine collisions:", [(v, s, sum(map(abs, v))) for v, s in gen])
for ell in d['examples'][:nex]:
    print("survivor", ell, "c=", [x % n for x in ell])
    for v, s in coll(ell):
        z = Zv(pis, v); D = Ds(z, s); nz = z[0] ** 2 + z[1] ** 2
        ch = 1
        for q in others: ch *= q ** vq(D, q)
        slack = 0.5 * math.log(nz * (2 if s % 2 else 1)) - math.log(ch)
        isgen = sum(a * b for a, b in zip(v, L)) % n == 0
        print(f"   v={v} norm={sum(map(abs, v))} s={s} genuine_lattice={isgen} D_s={D} factor={factor(abs(D)) if D else 0} "
              f"slack={slack:.2f} log q0={math.log(q0):.2f}")
