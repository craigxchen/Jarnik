"""Count, at each forced q0 of first4 (K given), the vectors v in S_{K-1} \\ Lambda_c (non-genuine at q0) and units s for
which ONE extra collision (v,s) at q0 is affordable (genuine labels on F \\ {q0}, all levels):
PU: q0 g_s <= nu_s^-1 sqrt(N);  ANG: J0 box meets circle with m_s = q0 g_s, others g;  C: q0 g_s <= |D_s|.
Also: does the genuine label map S_h (h = floor((K-1)/2)) onto Z/n (covering, rigtheory Thm 3.4)?
Usage: python3 affordable.py FFILE K p1..pr"""
import sys, math
from fractions import Fraction as Fr
from sllib import two_sq, genuine, Zv, Ds, VK, vq
from classify import load_F, j0_reason
Fpath, K = sys.argv[1], int(sys.argv[2]); ps = [int(x) for x in sys.argv[3:]]
F = load_F(Fpath); r = len(ps); pis = [two_sq(p) for p in ps]
h = (K-1)//2
Vd = VK(r, K-1); Sh = [tuple([0]*r)] + VK(r, h)
pre = {}
for v in Vd:
    z = Zv(pis, v); pre[v] = (z[0]**2+z[1]**2, [Ds(z, s) for s in range(4)])
print(f"ps={ps} K={K}: per forced q0: n, covered by S_{h}?, #nongenuine v (norm<=K-1, v~-v), affordable (v,s): PU / ANG / C")
QMAX = int(__import__("os").environ.get("QMAX", "10**9"))
for q0 in [q for q in F if q <= QMAX]:
    m, L = genuine(ps, q0); n = m//4
    cov = len(set(sum(a*b for a, b in zip(x, L)) % n for x in Sh)) == n
    cnt = [0, 0, 0]; nong = 0
    for v in Vd:
        if next(x for x in v if x) < 0: continue
        if sum(a*b for a, b in zip(v, L)) % n == 0: continue
        nong += 1
        N, D = pre[v]
        g = []
        for s in range(4):
            gg = 1
            for q in F:
                if q != q0: gg *= q**vq(D[s], q)
            g.append(gg)
        for s in range(4):
            ms = list(g); ms[s] *= q0
            pu = all(ms[t]**2 <= N*(2 if t % 2 else 1) for t in range(4))
            if not pu: continue
            cnt[0] += 1
            if j0_reason(ms, N) == 'ok': cnt[1] += 1
            if ms[s] <= abs(D[s]): cnt[2] += 1
    print(f"  q0={q0:4d} n={n:3d} cov={'Y' if cov else 'n'} nongen={nong:4d} affordable PU/ANG/C = {cnt[0]:4d} / {cnt[1]:4d} / {cnt[2]:4d}", flush=True)
