"""Anatomy of a full-reference survivor: extra collisions, D_s, |sin(arg Z_v - s pi/4)|."""
import sys, math
import deform, deform_full
from rlib14 import *
from forced import forced_set
K = int(sys.argv[1]); ps = [int(x) for x in sys.argv[2].split(",")]; q0 = int(sys.argv[3])
F = forced_set(ps, K, 60000); d = deform.run_q0(ps, K, F, q0, maxex=1)
pis = gaussian_primes(ps); m, n, L = d['m'], d['n'], d['L']
ell = d['examples'][0]; print("q0", q0, "survivor", ell, {k: d[k] for k in ['total', 'inclass', 'F1', 'F2']})
for v in VK_mod_sign(len(ps), K):
    t = sum(a * b for a, b in zip(v, ell)) % m
    if t % n: continue
    s = t // n; z = Zv(pis, v); D = Ds(z, s); ang = math.atan2(z[1], z[0]) - s * math.pi / 4
    gen = sum(a * b for a, b in zip(v, L)) % n == 0
    print(f"  v={v} s={s} genuine_lattice={gen} D_s={D} = {factor(abs(D))} |sin|={abs(math.sin(ang)):.2e} 1/q0={1/q0:.2e}")
