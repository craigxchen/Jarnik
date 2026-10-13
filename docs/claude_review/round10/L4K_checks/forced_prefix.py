# Tight set test: genuine labels (all levels) at forced q <= q1, then search good non-genuine labels at the
# forced q in (q1, q*].  Reports, for each q1, whether a good completion exists that is non-genuine (signature
# differs from genuine) at >= 1 / at all upper forced moduli, and its coherence on the upper forced moduli.
import sys, math, random, time
import numpy as np
import forced_coh
from forced_coh import *

def vq(x, q):
    x = abs(x); k = 0
    while x and x % q == 0: x //= q; k += 1
    return k

K, r, tag, q1, nsys, H = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6])
P, F = load_set(K, r, tag); V = monos(r, K)
B = np.array([[budget(P, v, s) for s in range(4)] for v in V])
pis = [gauss_prime(p) for p in P]
low = [q for q in F if q <= q1]; up = [q for q in F if q > q1]
load0 = np.zeros_like(B)
for i, v in enumerate(V):
    z = Zv(pis, v)
    for s in range(4):
        load0[i, s] = sum(vq(Ds(z, s), q) * math.log(q) for q in low)
assert (load0 <= B + 1e-9).all()
mods = [Mod(q, P, K, V) for q in up]
forced_coh.mods_by_q = {md.q: md for md in mods}
Zs = [(a, b) for a in range(1, math.isqrt(H) + 1) for b in range(0, math.isqrt(H) + 1) if a * a + b * b <= H]
Lz = {md.q: [md.T.L(z) if (z[0] ** 2 + z[1] ** 2) % md.q else -1 for z in Zs] for md in mods}
up7 = [q for q in up if q % 8 == 7]
rng = random.Random(5)
# monkeypatch: start load
orig = gen_good2
def gen_from(mods, B, rng):
    Bp = B - load0          # remaining budget
    return orig(mods, Bp, rng)
got = 0; fails = 0; coh_up = 0; coh_up7 = 0; nongen_all = 0; t0 = time.time()
while got < nsys and fails < 20:
    g = gen_from(mods, B, rng)
    if g is None: fails += 1; continue
    sysd, load = g; got += 1
    saved = globals().setdefault('saved', [])
    saved.append(dict(system={q: sysd[q][0] for q in sysd}, top_load=load))
    def isgen(md, ell):
        n = md.n; Lc = [x % n for x in md.L]
        return any(all((u * Lc[j] - ell[j]) % n == 0 for j in range(len(ell))) for u in range(1, n + 1) if math.gcd(u, n) == 1)
    ng = [md.q for md in mods if not isgen(md, sysd[md.q][0])]
    nongen_all += (len(ng) == len(mods))
    print(f"   system {got}: non-genuine-class at {len(ng)}/{len(mods)} upper moduli; non-genuine Q_7: {[q for q in ng if q % 8 == 7]}", flush=True)
    coh_up += coherent(sysd, mods, Zs, Lz, up); coh_up7 += coherent(sysd, mods, Zs, Lz, up7)
import pickle
pickle.dump(dict(P=P, F=F, K=K, q1=q1, load0=load0, saved=globals().get('saved', [])), open(f'prefix_{K}_{r}_{tag}_{q1}.pkl', 'wb'))
print(f"K={K} r={r} {tag} P={P}: genuine on forced q<={q1} ({low}); searched forced q in {up}: "
      f"good completions found {got} (failed searches {fails}); non-genuine class at every upper q: {nongen_all}; "
      f"coherent (norm<={H}) on all upper q: {coh_up}; on upper Q_7 {up7}: {coh_up7} [{time.time()-t0:.0f}s]", flush=True)
