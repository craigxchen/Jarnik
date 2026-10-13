"""Search for symmetric lattice polytopes A = {k in [-T,T]^d : |lam_j(k)| <= b_j} with rho(A) > zeta(A)
(the only place where a model can be faithful to the regime, compute.md Prop 2.2).
lam_j random small integer functionals (entries in {-1,0,1} by default, or {-2..2}).
rho mod p (upper bound for rho over Q: rank_p <= rank_Q); a hit is re-checked with exact rank."""
import sys, random, itertools
import flint
from core import monomials_lt, mono_val
from zono import zeta
from models import rho as rho_exact
from scan import rho_p

def build(d, T, lams):
    return [k for k in itertools.product(range(-T, T + 1), repeat=d)
            if all(abs(sum(a * b for a, b in zip(l, k))) <= bb for l, bb in lams)]

if __name__ == "__main__":
    d = int(sys.argv[1]); trials = int(sys.argv[2]); T = int(sys.argv[3]); cmax = int(sys.argv[4])
    seed = int(sys.argv[5]) if len(sys.argv) > 5 else 0
    rng = random.Random(seed)
    hits = 0
    for tr in range(trials):
        nl = rng.randint(1, 4)
        lams = []
        for _ in range(nl):
            l = tuple(rng.randint(-cmax, cmax) for _ in range(d))
            if not any(l):
                continue
            hw = sum(abs(x) for x in l) * T
            b = rng.randint(max(1, hw // 4), max(1, hw - 1))
            lams.append((l, b))
        A = build(d, T, lams)
        if len(A) < 2:
            continue
        r = rho_p(A)
        full = [(tuple(1 if y == x else 0 for y in range(d)), T) for x in range(d)] + lams
        z, Z = zeta(None, T, None, None, r, d=d, lams=full)
        tag = ""
        if z < r:
            re_ = rho_exact(A)
            tag = "GAP rho_exact=%d" % re_
            hits += 1
        print(tr, lams, "|A|=%d rho_p=%d zeta=%d %s" % (len(A), r, z, tag), flush=True)
    print("hits", hits)
