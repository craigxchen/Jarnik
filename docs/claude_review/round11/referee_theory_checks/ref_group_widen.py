"""Referee illustration (exact; far below the regime, evidence only): single-column widening versus the
referee's grouped widening on the round-5 Paley fake profile at q = 7 (n = 8, r = 64).
Box: T = 1, narrow slabs |<k,beta_q>| <= w = 2 for q not widened; widened q: |l_q| <= lam*m/2.
UM_p tested exactly for every k* (sl_exact.py), p = 3.
Prediction (proved in the report): single widening of an unflipped column p fails, because its 7 duplicate
copies keep |l_p| <= w/2, so l_p has <= w+1 = 3 < m levels."""
import itertools, sys
from sl_exact import um_points, certificate
from ref_fake_sl import round5_profile

def box_model(cols, n, T, w, m, widened, lam):
    d = n - 1
    S = []; L = {q: [] for q in widened}
    keyset = {}
    # deduplicate identical narrow constraints for speed
    narrow = []
    seen = set()
    for qi, col in enumerate(cols):
        if qi in widened:
            continue
        key = tuple(col[x] - col[d] for x in range(d))
        if key not in seen:
            seen.add(key); narrow.append(key)
    wid = {q: tuple(cols[q][x] - cols[q][d] for x in range(d)) for q in widened}
    for k in itertools.product(range(-T, T + 1), repeat=d):
        if any(abs(sum(a * b for a, b in zip(k, key))) > w for key in narrow):
            continue
        ok = True
        lv = {}
        for q, key in wid.items():
            val = sum(a * b for a, b in zip(k, key))
            assert val % 2 == 0
            if abs(val // 2) > lam * m // 2:
                ok = False; break
            lv[q] = val // 2
        if ok:
            S.append(k)
            for q in widened:
                L[q].append(lv[q])
    return S, L

if __name__ == '__main__':
    cols, n = round5_profile(7)
    b = 8
    # columns 0..7 are the 8 unflipped copies of label 1; flipped copies follow after b*(n-1) columns
    lab1 = list(range(0, b))
    flipped_lab1 = [b * (n - 1) + x for x in range(n) if (x if x != 0 else 1) == 1]   # rows with a(x) = 1
    T, w = 1, 2
    for m in (4, 5):
        for desc, widened, lam in (("single unflipped column 0", [0], 1),
                                   ("group: label-1 class (8 copies + flipped)", lab1 + flipped_lab1, 1),
                                   ("group, lam=2", lab1 + flipped_lab1, 2)):
            S, L = box_model(cols, n, T, w, m, widened, lam)
            p = 3
            lv = L[0]
            dimV, good = um_points(S, m, lv, p) if len(S) > 0 else (0, [])
            levs = sorted(set(lv)) if lv else []
            ver = ""
            if good:
                ok, c = certificate(S, m, lv, p, good[0]); ver = f" certificate verified: {ok}"
            bylev = {}
            for i in good:
                bylev[lv[i]] = bylev.get(lv[i], 0) + 1
            print(f"m={m} {desc}: |S|={len(S)} levels of l_p: {levs} dimV={dimV} UM_p points by level {bylev}{ver}")
            sys.stdout.flush()
