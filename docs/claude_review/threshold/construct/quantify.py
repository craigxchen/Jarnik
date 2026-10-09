"""Table for construct.md section 2.7: for named measures c, the order, the cube-profile payoff
MW(supp c) = E_eps h_S(eps) (eps uniform in {+-1}^M), the uniform-balanced payoff 2 E_bal w, and the best
half-separating adversary 2*adv (exact LP).  The method of route (C) needs ord > every payoff."""
import itertools, random
from fractions import Fraction as Fr
from clib import widths, adversary, popc
from gens import orbit, Qslice
from peeldp import bound

def payoffs(S, M):
    w = widths(S, M)
    mw = Fr(2 * 2 * sum(w.values()), 2**M)                 # 2 E_iid w (complement pairs counted once -> x2)
    q = M // 2
    bal = [v for P, v in w.items() if min(popc(P), M - popc(P)) == q]
    eb = 2 * Fr(sum(bal), len(bal))
    adv, _ = adversary(w, M)
    return mw, eb, 2 * adv

rows = []
for M in range(3, 8):
    r = (M - 1) // 2
    v = tuple(range(-r, M - r))
    S = orbit(v)
    ordc = M * (M - 1) // 2
    rows.append((f"Vandermonde alternant M={M}, v={v}", M, ordc) + payoffs(S, M))
# generalized alternant on S^5(6,3), order 14 (genalt.py)
reps = [(5,1,0,-1,-2), (4,2,0,-1,-2), (3,2,1,0,-3), (3,2,1,-1,-2)]
S = sorted(set(itertools.chain.from_iterable(orbit(x) for x in reps)))
rows.append(("generalized alternant on S^5(6,3)", 5, 14) + payoffs(S, 5))
# Q^M(p,p) slice: reg = floor(phi) exact for these (upper.md); take the alternant-attaining points
for M, p, regv in [(4, 2, 6), (5, 3, 10), (6, 2, 6)]:
    S = Qslice(M, p, p)
    rows.append((f"full slice Q^{M}({p},{p}) (reg={regv})", M, regv) + payoffs(S, M))
for name, M, o, mw, eb, adv2 in rows:
    print(f"{name:45s} ord={o:3d}  MW={str(mw):7s} ord/MW={float(o/mw):.4f}  2E_bal w={str(eb):7s}  2adv={str(adv2):7s} ord/(2adv)={float(o/adv2):.4f}")

# spike profile M=7 (symmetric cut polytope), rigorous peeling upper bound for reg
from clib import sym_adversary
M = 7
for n in (2, 4, 8):
    U = tuple(n * u for u in (6, 5, 4, 4, 5, 6))
    pb = bound(M, U, 0)
    w = {1: 12 * n, 2: 10 * n, 3: 8 * n}
    from math import comb
    mw = Fr(2 * sum(comb(M, j) * w[min(j, M - j)] for j in range(1, M)), 2**M)
    adv, _ = sym_adversary(w, M)
    print(f"spike profile M=7, n={n}: reg <= {pb} (peeling)  MW={mw} ({float(pb/mw):.4f})  2adv={2*adv} ({float(pb/(2*adv)):.4f})")

# per-prime inequality of Proposition C4:  T(b)+T(-b) <= beta_u + beta_v - e*h_S(eps)
random.seed(3)
bad = 0
for t in range(20000):
    M = random.randint(2, 5); e = random.randint(1, 4)
    S = list({tuple(random.randint(-3, 3) for _ in range(M)) for _ in range(random.randint(1, 8))})
    al = {k: random.choice([0, 0, 1, 2, 5]) for k in S}; alb = {k: random.choice([0, 0, 1, 3]) for k in S}
    eps = [random.choice([1, -1]) for _ in range(M)]
    def x(k, sg): return Fr(e, 2) * sg * sum(a * b for a, b in zip(k, eps))
    T = lambda sg: min(al[k] + x(k, sg) for k in S) + min(alb[k] - x(k, sg) for k in S)
    u = min(S, key=lambda k: x(k, 1)); v = max(S, key=lambda k: x(k, 1))
    h = x(v, 1) - x(u, 1)
    if T(1) + T(-1) > al[u] + alb[u] + al[v] + alb[v] - 2 * h: bad += 1
print("Proposition C4 per-prime inequality: violations in 20000 random instances:", bad)
