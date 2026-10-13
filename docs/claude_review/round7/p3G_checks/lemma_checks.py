"""Brute-force checks of the elementary lemmas of p3G.md section 2 and 5 (floating point; evidence only).
(a) Lemma 2.2: for tight (v,s) the admissible set of x = <v,phi-theta> mod pi is sigma*[0, pi-2d].
(b) Lemma 2.3 (length-2 and self circuits): admissible x is {0} when the arc condition holds.
(c) Lemma 5.2 (two angles in one half, units 3 and 1): opposite orientations and d + d' > pi/2.
"""
import math, random
rng = random.Random(7)
pi = math.pi
def dist(t):  # distance to pi Z
    t = t % pi
    return min(t, pi - t)
def orient(alpha, s):
    a = (alpha - s * pi / 4) % pi
    return (1 if a < pi / 2 else -1), min(a, pi - a)
bad = 0
for _ in range(20000):
    alpha = rng.uniform(0, 2 * pi); s = rng.randrange(4)
    sg, d = orient(alpha, s)
    for x in [rng.uniform(-pi, pi) for _ in range(20)]:
        adm = dist(alpha + x - s * pi / 4) >= d - 1e-12
        xm = (sg * x) % pi
        pred = xm <= pi - 2 * d + 1e-12 or xm >= pi - 1e-12
        if adm != pred and abs(dist(alpha + x - s * pi / 4) - d) > 1e-9:
            bad += 1
print('(a) mismatches:', bad)
bad = 0; tested = 0
for _ in range(20000):
    a1 = rng.uniform(0, 2 * pi); a2 = rng.uniform(0, 2 * pi); s1 = rng.randrange(4); s2 = rng.randrange(4)
    sg1, d1 = orient(a1, s1); sg2, d2 = orient(a2, s2)
    if sg1 == sg2 or d1 + d2 <= pi / 2:
        continue
    tested += 1
    # same monomial direction: x1 = x2 = x (circuit sigma1 v + sigma2 v = 0)
    for x in [rng.uniform(-pi, pi) for _ in range(50)]:
        ok1 = dist(a1 + x - s1 * pi / 4) >= d1 - 1e-12
        ok2 = dist(a2 + x - s2 * pi / 4) >= d2 - 1e-12
        if ok1 and ok2 and dist(x) > 1e-9:
            bad += 1
print('(b) pair circuits tested', tested, 'admissible nonzero x found:', bad)
bad = 0
for _ in range(100000):
    h = rng.randrange(2)
    lo = pi / 4 if h == 0 else 3 * pi / 4
    al, ar = sorted([rng.uniform(lo, lo + pi / 2) for _ in range(2)])
    if h == 0:
        (sa, da), (sb, db) = orient(al, 3), orient(ar, 1)
    else:
        (sa, da), (sb, db) = orient(al, 1), orient(ar, 3)
    if not (sa == -sb and da + db > pi / 2):
        bad += 1
print('(c) failures:', bad)
