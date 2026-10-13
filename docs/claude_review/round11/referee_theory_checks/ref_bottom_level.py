"""Referee check: on supports where the level function l takes values in m+1 (or fewer) consecutive integers,
UM_p can only occur at the bottom level (exactly m+1 levels) and never with <= m levels.
Exact (sl_exact.py).  Random small models in dimension d = 2, 3 with random extra slabs.
"""
import random, itertools, sys
from sl_exact import um_points

def model(d, T, m, nslab, wslab, rng, lvec, half):
    # S = box [-T,T]^d cut by |<k,u_q>| <= wslab (random u_q) and the l-slab l(k) in [-half, half]
    us = [[rng.choice((-1, 0, 1)) for _ in range(d)] for _ in range(nslab)]
    S = []
    for k in itertools.product(range(-T, T + 1), repeat=d):
        l = sum(a * b for a, b in zip(k, lvec))
        if abs(l) > half:
            continue
        if any(abs(sum(a * b for a, b in zip(k, u))) > wslab for u in us):
            continue
        S.append(k)
    return S

if __name__ == '__main__':
    rng = random.Random(20261013)
    tot = {'odd_um': 0, 'even_nonbottom_um': 0, 'even_bottom_um': 0, 'instances': 0, 'even_inst_with_um': 0}
    for trial in range(60):
        d = rng.choice((2, 3))
        T = rng.choice((2, 3)) if d == 3 else rng.choice((3, 4, 5))
        m = rng.choice((2, 3, 4, 5, 6)) if d == 2 else rng.choice((2, 3, 4))
        lvec = [rng.choice((1, 2, -1, 1)) for _ in range(d)]
        half = m // 2                     # |l| <= m/2 : the theory.md widening
        S = model(d, T, m, rng.choice((0, 1, 2)), rng.choice((1, 2, 3)), rng, lvec, half)
        if len(S) < 3:
            continue
        lv = [sum(a * b for a, b in zip(k, lvec)) for k in S]
        p = rng.choice((2, 3, 5))
        dimV, good = um_points(S, m, lv, p)
        tot['instances'] += 1
        levels = sorted(set(lv))
        bottom = min(lv)
        nb = [i for i in good if lv[i] != bottom]
        if m % 2 == 1:
            tot['odd_um'] += len(good)
        else:
            tot['even_nonbottom_um'] += len(nb)
            tot['even_bottom_um'] += len(good) - len(nb)
            if good:
                tot['even_inst_with_um'] += 1
        print(f"d={d} T={T} m={m} p={p} |S|={len(S)} levels={levels[0]}..{levels[-1]} ({len(levels)}) dimV={dimV} "
              f"UM points={len(good)} (non-bottom {len(nb)})")
        sys.stdout.flush()
    print("TOTALS:", tot)
