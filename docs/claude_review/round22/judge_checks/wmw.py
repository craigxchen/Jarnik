"""Judge check J3 (evidence only): test the weighted mean-width inequality proposed in toric.md open point 1,
     reg(S) <= (1/6) * sum_{|A|=2} width_S(1_A)            (W)
for finite S in M = {x in Z^5 : sum x = 0}. reg(S) = least d such that polynomials of degree <= d restricted to S
span Q^S (= max ord_e of a Laurent polynomial with support in S; toric.md Lemma 3.1). Ranks are computed mod p;
full rank mod p implies full rank over Q, so the computed reg is an UPPER bound for the rational reg
(rank mod p <= rank over Q). Hence a violation of (W) found here must be re-certified over Q; a non-violation
is only evidence.
Sanity: also checks construct.md Thm B, reg(S) <= 2^{-4} sum_{A subset [5]} width_S(1_A) (refereed), which must hold.
"""
import itertools, random, sys
from fractions import Fraction
P = 2147483647

def rank_mod(rows):
    rows = [r[:] for r in rows]
    rk = 0; ncol = len(rows[0]) if rows else 0
    for c in range(ncol):
        piv = None
        for i in range(rk, len(rows)):
            if rows[i][c] % P:
                piv = i; break
        if piv is None: continue
        rows[rk], rows[piv] = rows[piv], rows[rk]
        inv = pow(rows[rk][c], P - 2, P)
        rows[rk] = [(v * inv) % P for v in rows[rk]]
        for i in range(len(rows)):
            if i != rk and rows[i][c] % P:
                f = rows[i][c]
                rows[i] = [(a - f * b) % P for a, b in zip(rows[i], rows[rk])]
        rk += 1
    return rk

def monos(d, n=4):
    out = []
    for e in itertools.product(range(d + 1), repeat=n):
        if sum(e) <= d: out.append(e)
    return out

def reg(S):
    pts = [s[:4] for s in S]
    n = len(pts)
    d = 0
    while True:
        cols = monos(d)
        rows = [[(pow(p[0], e[0]) * pow(p[1], e[1]) * pow(p[2], e[2]) * pow(p[3], e[3])) % P for e in cols] for p in pts]
        if rank_mod(rows) == n: return d
        d += 1

subsets2 = list(itertools.combinations(range(5), 2))
allA = [A for r in range(0, 6) for A in itertools.combinations(range(5), r)]

def width(S, A):
    vals = [sum(s[a] for a in A) for s in S]
    return max(vals) - min(vals)

def wbound(S):
    return Fraction(sum(width(S, A) for A in subsets2), 6)

def thmB(S):
    return Fraction(sum(width(S, A) for A in allA), 16)

def rand_pt(c):
    while True:
        x = [random.randint(-c, c) for _ in range(4)]
        x5 = -sum(x)
        if abs(x5) <= c: return tuple(x + [x5])

def main(seed, trials):
    random.seed(seed)
    worst = Fraction(0); viol = 0; tight = 0; bviol = 0
    for t in range(trials):
        c = random.choice([1, 1, 2, 2, 3])
        n = random.randint(2, 26)
        S = list({rand_pt(c) for _ in range(n)})
        if len(S) < 2: continue
        r = reg(S); b = wbound(S); tb = thmB(S)
        if r > tb: bviol += 1
        if r > b:
            viol += 1
            print("VIOLATION", r, b, S); sys.stdout.flush()
        if r == b: tight += 1
        if Fraction(r) / b > worst: worst = Fraction(r) / b
    print(f"seed {seed}: trials {trials}, (W) violations {viol}, ThmB violations {bviol}, tight {tight}, max reg/bound {worst} = {float(worst):.4f}")

if __name__ == "__main__":
    main(int(sys.argv[1]), int(sys.argv[2]))
