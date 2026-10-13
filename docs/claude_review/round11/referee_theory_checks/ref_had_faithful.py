"""Referee check (exact): SL-full on the Sylvester-Hadamard profile n = 8 with the widening theory.md actually
uses, |<k,beta_p>| <= m/e_p (i.e. |l_p| <= m/2), versus the evidence box of theory_checks/had_sl.py,
|<k,beta_p>| <= 2m (|l_p| <= m).  T = 1, narrow slabs |<k,beta_q>| <= 2 (q != p).  All k* tested (no sampling),
p = 3 and p = 5.  Exact rational nullspace + p-saturation (sl_exact.py)."""
import itertools, sys
from sl_exact import um_points, certificate

def build(m, widen):
    pts = list(itertools.product((0, 1), repeat=3))
    n = len(pts); d = n - 1
    cols = [[(-1) ** (sum(a * b for a, b in zip(q, x)) % 2) for x in pts] for q in pts[1:]]
    p_idx = 0
    S = []; lv = []
    for k in itertools.product((-1, 0, 1), repeat=d):
        ok = True
        lp = None
        for qi, col in enumerate(cols):
            val = sum(k[x] * (col[x] - col[d]) for x in range(d))   # <k,beta_q> on the zero-sum slice
            if qi == p_idx:
                if abs(val) > widen:
                    ok = False; break
                lp = val // 2                                     # beta = 2a - 1 => <k,beta> = 2 l
                assert val % 2 == 0
            elif abs(val) > 2:
                ok = False; break
        if ok:
            S.append(k); lv.append(lp)
    return S, lv

if __name__ == '__main__':
    for m in (3, 4, 5):
        for label, widen in (("theory A_p: |<k,b_p>|<=m", m), ("had_sl box: |<k,b_p>|<=2m", 2 * m)):
            S, lv = build(m, widen)
            for p in (3, 5):
                dimV, good = um_points(S, m, lv, p)
                levs = sorted(set(lv))
                bylev = {}
                for i in good:
                    bylev[lv[i]] = bylev.get(lv[i], 0) + 1
                cnt = {l: lv.count(l) for l in levs}
                ver = ""
                if good:
                    ok, c = certificate(S, m, lv, p, good[0])
                    ver = f"; certificate for first UM point verified: {ok}"
                print(f"m={m} {label} p={p}: |S|={len(S)} levels {levs[0]}..{levs[-1]} sizes {cnt} dimV={dimV} "
                      f"UM points by level {bylev}{ver}")
                sys.stdout.flush()
