# Referee check of the class/digit mechanism of X.2 and of the word counts.
# Own implementation; exact integers only.
import random, itertools
from math import comb, ceil, log2, gcd
from functools import lru_cache
random.seed(1)

def classsum_vectors(b, H):
    """all nonzero integer vectors in Z^b with l1 mass <= H"""
    out = []
    def rec(i, rem, cur):
        if i == b:
            if any(cur): out.append(tuple(cur))
            return
        for x in range(-rem, rem+1):
            cur.append(x); rec(i+1, rem-abs(x), cur); cur.pop()
    rec(0, H, [])
    return out

def detecting_rows(b, H):
    T = ceil(2*H*log2(2*b)) + 2
    vecs = classsum_vectors(b, H)
    while True:
        rows = [[random.randint(0, 1) for _ in range(b)] for _ in range(T)]
        if all(any(sum(c*x for c, x in zip(r, v)) != 0 for r in rows) for v in vecs):
            return T, rows, len(vecs)

def check_cells(H, b, s, D_exp_extra=3, trials=3000):
    B = 64*H*H
    T, rows, nv = detecting_rows(b, H)
    BT = B**T
    g = 2*B; D = 3
    while D < 100*g*BT: D *= 3   # D = r^{k_r} with r = 3, D/B^T >= 100 g
    cls = [j % b for j in range(s)]
    def cell_value(j):
        z = 0
        for t in range(T):
            a = (B//(4*H)) * rows[t][cls[j]]
            z = z*B + a + random.randint(0, 1)
        lo = (D*z + BT - 1)//BT; hi = (D*(z+1))//BT
        return random.randint(lo, hi)
    L = [cell_value(j) for j in range(s)]
    # (i) digit window: B^(t-1) L_j / D mod 1 lies in [a_t/B, (a_t+2)/B]
    for j in range(s):
        for t in range(1, T+1):
            a = (B//(4*H))*rows[t-1][cls[j]]
            fr = (B**(t-1)*L[j]) % D          # fractional part times D
            assert a*D <= fr*B <= (a+2)*D, (j, t)
    # (ii) tested words: D does not divide sum v_j L_j
    tested = failed = 0
    for _ in range(trials):
        h = random.randint(1, H)
        v = {}
        for _ in range(h):
            j = random.randrange(s); v[j] = v.get(j, 0) + random.choice((1, -1))
        v = {j: c for j, c in v.items() if c}
        if not v: continue
        cs = [0]*b
        for j, c in v.items(): cs[cls[j]] += c
        if not any(cs): continue
        tested += 1
        S = sum(c*L[j] for j, c in v.items())
        if S % D == 0: failed += 1
        # quantitative: distance of B^(t-1) S / D from Z is >= 1/(4H) - 2H/B for detecting t
        t = next(t for t in range(T) if sum(rows[t][k]*cs[k] for k in range(b)) != 0)
        fr = (B**t*S) % D
        dist = min(fr, D-fr)
        assert dist*B*4*H >= D*(B - 8*H*H), "distance bound"
    # (iii) same-class difference words: no control; show they can be divisible
    return T, nv, tested, failed

for (H, b, s) in [(2, 3, 60), (2, 4, 80), (3, 3, 60)]:
    T, nv, tested, failed = check_cells(H, b, s)
    print(f"H={H} b={b}: T={T} rows detect all {nv} class-sum vectors; "
          f"{tested} random tested words, D | sum in {failed} cases")

# (iv) exact counts of words of mass h with all class sums zero
@lru_cache(None)
def Z(n, h, total=0):
    """number of v in Z^n, ||v||_1 = h, sum v = total"""
    if n == 0: return 1 if (h == 0 and total == 0) else 0
    return sum(Z(n-1, h-abs(x), total-x) for x in range(-h, h+1))
def allwords(s, h):
    return sum(comb(s, k)*comb(h-1, k-1)*2**k for k in range(1, h+1)) if h else 1
def excluded(sizes, h):
    poly = [1] + [0]*h
    for n in sizes:
        zn = [Z(n, k) for k in range(h+1)]
        poly = [sum(poly[i]*zn[k-i] for i in range(k+1)) for k in range(h+1)]
    return poly[h]
print("\nfraction of signed words of mass h with all class sums zero (balanced classes):")
for s, b in [(240, 4), (240, 8), (240, 16), (960, 32)]:
    sizes = [s//b + (1 if i < s % b else 0) for i in range(b)]
    row = []
    for h in range(1, 7):
        e, a = excluded(sizes, h), allwords(s, h)
        row.append(f"h={h}:{e/a:.5f}")
    diff_frac = sum(comb(n, 2) for n in sizes)/comb(s, 2)
    print(f" s={s} b={b}: " + " ".join(row) + f" | same-class share of differences {diff_frac:.5f} (1/b={1/b:.5f})")
# odd mass: all words tested (total coefficient sum = mass mod 2)
# Turan: a test graph on differences avoiding every intra-class pair is b-partite,
# hence has at most (1-1/b) s^2/2 edges.
