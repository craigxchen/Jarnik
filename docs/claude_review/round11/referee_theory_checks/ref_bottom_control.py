"""Controls for sl_exact.py: (1) with a WIDER l-range (|l| <= m, 2m+1 levels) non-bottom UM points do occur,
so the bottom-level phenomenon is not an artefact of the tester; (2) every reported UM point (in a sample)
is confirmed by an explicit certificate c, verified from scratch (M c = 0 exactly; unique p-adic minimiser)."""
import random, itertools
from sl_exact import um_points, certificate
from ref_bottom_level import model

rng = random.Random(7)
nonbottom = 0; certs_ok = 0; certs = 0
for trial in range(25):
    d = 2; T = rng.choice((3, 4)); m = rng.choice((3, 4, 5)); lvec = [rng.choice((1, 2)), rng.choice((1, -1))]
    S = model(d, T, m, rng.choice((0, 1)), rng.choice((2, 3)), rng, lvec, m)   # |l| <= m: wide
    if len(S) < 3:
        continue
    lv = [lvec[0] * k[0] + lvec[1] * k[1] for k in S]
    p = rng.choice((2, 3))
    dimV, good = um_points(S, m, lv, p)
    nb = [i for i in good if lv[i] != min(lv)]
    nonbottom += len(nb)
    for i in good[:3]:
        ok, c = certificate(S, m, lv, p, i)
        certs += 1; certs_ok += ok
    print(f"wide: T={T} m={m} p={p} |S|={len(S)} levels {min(lv)}..{max(lv)} dimV={dimV} UM={len(good)} non-bottom={len(nb)}")
# certificates in the faithful (m+1 levels) regime too
rng = random.Random(20261013)
for trial in range(15):
    d = 2; T = rng.choice((3, 4)); m = rng.choice((4, 6)); lvec = [rng.choice((1, 2)), rng.choice((1, -1))]
    S = model(d, T, m, rng.choice((0, 1)), rng.choice((2, 3)), rng, lvec, m // 2)
    if len(S) < 3:
        continue
    lv = [lvec[0] * k[0] + lvec[1] * k[1] for k in S]
    p = rng.choice((2, 3))
    dimV, good = um_points(S, m, lv, p)
    for i in good[:3]:
        ok, c = certificate(S, m, lv, p, i)
        certs += 1; certs_ok += ok
print(f"non-bottom UM points in wide models: {nonbottom}; explicit certificates verified: {certs_ok}/{certs}")
