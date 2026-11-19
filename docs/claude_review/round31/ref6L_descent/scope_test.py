"""Astra's scope example G_e = B + c_e, plus random affine pencils G_e = a_e B + c_e.
Reports H, log Delta / log H, and nu = 6 + log Delta / log P (P = H^(1/20)).
Also checks exactness of t_ij (integrality after division by |G_ij|^2)."""
import random, math, sys
sys.path.insert(0, '.')
from gauss import *

random.seed(7)

def report(G, label):
    z = tuple_from_blocks(G)
    H, logD = width_and_H(z)
    logH = math.log(H)
    t = cofactors(G)
    tmax = max(abs(v) for v in t.values())
    print(f"{label}: log10 H={logH/math.log(10):.1f}  logD/logH={logD/logH:+.4f}  nu={6+20*logD/logH:+.3f}  max|t|={tmax}  min|t|={min(abs(v) for v in t.values())}")

# scope example
for trial in range(3):
    while True:
        c = [(random.randint(-5, 5), random.choice([-4,-3,-2,-1,1,2,3,4])) for _ in range(10)]
        cc = set(c) | set(conj(x) for x in c)
        if len(cc) == 20:
            break
    for B in [10**3, 10**6, 10**12, 10**24]:
        G = [(B+x[0], x[1]) for x in c]
        report(G, f"scope c#{trial} B=1e{int(math.log10(B))}")

# random pencils with aligned a_e = Gaussian units times small real (keep phases aligned) and Gaussian B
for trial in range(3):
    a = [(random.randint(1, 3), 0) for _ in range(10)]
    c = [(random.randint(-6, 6), random.randint(-6, 6)) for _ in range(10)]
    for s in [10**4, 10**8, 10**16]:
        B = (s + random.randint(0, s), random.randint(0, s))
        G = [(ai[0]*B[0]+ci[0], ai[0]*B[1]+ci[1]) for ai, ci in zip(a, c)]
        report(G, f"pencil#{trial} |B|~{s:.0e}")
# unaligned a_e
for trial in range(2):
    a = [(random.randint(-3, 3), random.randint(-3, 3)) for _ in range(10)]
    a = [x if x != (0, 0) else (1, 0) for x in a]
    c = [(random.randint(-6, 6), random.randint(-6, 6)) for _ in range(10)]
    B = (10**8+random.randint(0, 10**8), random.randint(0, 10**8))
    G = [mul(ai, B) for ai in a]
    G = [(g[0]+ci[0], g[1]+ci[1]) for g, ci in zip(G, c)]
    report(G, f"unaligned#{trial}")
