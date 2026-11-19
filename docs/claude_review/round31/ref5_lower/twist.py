"""Unit twists at deg T < d (WLOG T = 0, Y = theta): raw F_S has leading coefficient theta or theta^2
(the side contains exactly one of T, T+1).  The record's model uses MONIC points, so the natural points
are F_S / lc(F_S).  Compare pair-difference degrees untwisted vs monic-normalised, and search clusters.
Also: pairs of opposite sides of one form (z and zbar), excluded from Astra's statement."""
import sys, itertools
from f4poly import *
from enum_gcd import prime_data

d = int(sys.argv[1])
T = int(sys.argv[2]) if len(sys.argv) > 2 else 0
Q = 1 << d
fs = forms_and_sides(d)
sides = []
for fi, (w, k, c) in enumerate(fs):
    for S in (k, c):
        F = prod(linfac(T, v) for v in sorted(S))
        sides.append((fi, S, F))
D0 = deg(sides[0][2])
print(f"d={d} T={T}: D0={D0}, number of sides {len(sides)}")
lcs = {}
for fi, S, F in sides:
    lcs[F[-1]] = lcs.get(F[-1], 0) + 1
print("leading coefficients of raw F_S:", lcs)
mon = [trim([MUL[INV[F[-1]]][x] for x in F]) for (_, _, F) in sides]

def stats(points, label):
    best = {}
    for x, y in itertools.combinations(range(len(sides)), 2):
        fx, fy = sides[x][0], sides[y][0]
        dd = deg(add(points[x], points[y]))
        kind = "same form (z,zbar)" if fx == fy else "distinct forms"
        best.setdefault(kind, []).append(dd)
    for kind, L in best.items():
        print(f"  {label:9s} {kind:20s}: min deg diff {min(L)}, max {max(L)}, D0/2 = {D0/2}")

stats([F for (_, _, F) in sides], "raw")
stats(mon, "monic")

# cluster search on monic points (any set of sides, forms may repeat), full common gcd division
nS = len(sides)
dif = [[deg(add(mon[x], mon[y])) if x != y else -1 for y in range(nS)] for x in range(nS)]
best = {}
def visit(i, chosen, G, r):
    if len(chosen) >= 2:
        g = deg(G)
        if r <= g:
            s = len(chosen)
            pd = D0 - g
            if s not in best or pd < best[s][0]:
                best[s] = (pd, tuple(chosen))
        else:
            return
    for j in range(i, nS):
        nG = gcd(G, mon[j]) if chosen else mon[j]
        nr = max([r] + [2 * dif[x][j] - D0 for x in chosen])
        if chosen and nr > deg(nG):
            continue
        visit(j + 1, chosen + [j], nG, nr)
visit(0, [], (), -10**9)
print("monic-normalised clusters at this T (any sides, full common gcd):")
for s in sorted(best):
    print(f"  s={s}: min point degree {best[s][0]}  sides {best[s][1]}")
