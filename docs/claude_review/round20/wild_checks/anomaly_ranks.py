# Exact rank computations for the anomaly lemma (Thm A of wild.md).
# (1) Character lattice of M_{0,k} in Pluecker coordinates: rank of all cross-ratio exponent vectors
#     [ab;cd] = p_ac p_bd / (p_ad p_bc) in Z^{binom(k,2)}; expected k(k-3)/2 (= dim of Pluecker torus mod column scaling).
# (2) Rank of the cross-ratios supported on the first s points (s=4,5): expected s(s-3)/2.
# (3) Anomaly inequality for forgetful fibres Y (fibre of M_{0,k} -> M_{0,s}):
#     dim Y = k-s, coset codim = s(s-3)/2, D = k(k-3)/2, dim V = k-3; anomalous iff dimY>=1 and dimY > max(0, dimV + (D - codim) - D).
# (4) Chord linear space L_M: dim V = M-2 in torus of dim binom(M,2)-1; triangle-shape fibres dim M-3, coset codim 2.
from itertools import combinations
from fractions import Fraction

def rank(rows):
    rows = [list(map(Fraction, r)) for r in rows]
    if not rows: return 0
    m = len(rows[0]); r = 0
    for c in range(m):
        piv = None
        for i in range(r, len(rows)):
            if rows[i][c] != 0: piv = i; break
        if piv is None: continue
        rows[r], rows[piv] = rows[piv], rows[r]
        for i in range(len(rows)):
            if i != r and rows[i][c] != 0:
                f = rows[i][c] / rows[r][c]
                rows[i] = [a - f*b for a, b in zip(rows[i], rows[r])]
        r += 1
    return r

def crossratio_vectors(k, pts):
    pairs = list(combinations(range(k), 2)); idx = {p: i for i, p in enumerate(pairs)}
    def e(a, b): return idx[(min(a, b), max(a, b))]
    vecs = []
    for a, b, c, d in combinations(pts, 4):
        for (A, B, Cc, Dd) in [(a, b, c, d), (a, c, b, d), (a, d, b, c)]:
            v = [0]*len(pairs)
            v[e(A, Cc)] += 1; v[e(B, Dd)] += 1; v[e(A, Dd)] -= 1; v[e(B, Cc)] -= 1
            vecs.append(v)
    return vecs

if __name__ == "__main__":
    for k in range(4, 10):
        full = rank(crossratio_vectors(k, range(k)))
        D = k*(k-3)//2
        line = f"k={k}: rank(all cross-ratios)={full} (expected {D})"
        for s in (4, 5):
            if s <= k:
                rs = rank(crossratio_vectors(k, range(s)))
                dimY = k - s; dimV = k - 3
                thr = max(0, dimV + (D - rs) - D)
                anom = dimY >= 1 and dimY > thr
                line += f" | s={s}: rank={rs} (exp {s*(s-3)//2}), dimY={dimY}, threshold={thr}, anomalous={anom}"
        print(line)

    print()
    for M in range(3, 12):
        n = M*(M-1)//2 - 1; dimV = M - 2; dimY = M - 3; codim = 2
        thr = max(0, dimV + (n - codim) - n)
        print(f"chord space M={M}: dimV={dimV}, torus dim={n}, triangle fibre dimY={dimY}, threshold={thr}, anomalous={dimY>=1 and dimY>thr}")
