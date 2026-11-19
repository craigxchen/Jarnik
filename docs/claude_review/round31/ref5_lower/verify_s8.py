"""Explicit check of the d=4, s=8 cluster of point degree 104 found by family_fast.py (deg T = 15):
build the 8 marks, divide by the Euclidean gcd, check norms and every pair difference."""
from f4poly import *

d = 4
Q = 16
fs = forms_and_sides(d)
sides_idx = (0, 2, 5, 7, 8, 10, 13, 15)
found = None
for low in range(1 << (15 - d)):
    T = (1 << 15) | (low << d)
    Fs = []
    for sidx in sides_idx:
        fi, ch = divmod(sidx, 2)
        S = fs[fi][1 + ch]
        Fs.append(prod(linfac(T, v) for v in sorted(S)))
    G = Fs[0]
    for F in Fs[1:]:
        G = gcd(G, F)
    pts = [divmod_(F, G)[0] for F in Fs]
    D = deg(pts[0])
    if D != 104:
        continue
    N = mul(pts[0], conj(pts[0]))
    assert all(x in (0, 1) for x in N)
    assert all(mul(p, conj(p)) == N for p in pts)
    mx = max(deg(add(pts[i], pts[j])) for i in range(8) for j in range(i))
    if 2 * mx <= D:
        found = (T, D, mx, deg(G), deg(N))
        break
print("T =", bin(found[0]), "point degree", found[1], "max pair difference degree", found[2],
      "gcd degree", found[3], "parent degree", found[4])
print("Astra V.1 at m=8: point degree 112, parent 224")
