# EXACT (mod p ranks). Test of the cross-polytope criticality conjecture  reg(S) <= Omega'(S) := E_eps h_S(<eps,.>)
# (eps uniform in {+-1}^rho) on random finite S in Z^rho, and comparison with Theorem B's MW(S) (threshold/construct.md).
# reg(S) = least d such that every function on S is the restriction of a polynomial of degree <= d.
import itertools, random, sys, math, numpy as np
from common import monomials, rank_rows_mod_p, P31
random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
def reg(S):
    rho = len(S[0]); d = 0
    while True:
        rows = np.array([[math.prod(s[i]**a[i] for i in range(rho)) % P31 for s in S] for a in monomials(rho, d+1)], dtype=np.int64)
        if rank_rows_mod_p(rows) == len(S): return d
        d += 1
def width(S, f): v = [sum(a*b for a, b in zip(f, s)) for s in S]; return max(v)-min(v)
def omega_cross(S):
    rho = len(S[0]); return sum(width(S, e) for e in itertools.product([1, -1], repeat=rho)) / 2**rho
def MW(S):
    rho = len(S[0]); return 2*sum(width(S, [1 if i in P else 0 for i in range(rho)]) for k in range(rho+1) for P in itertools.combinations(range(rho), k)) / 2**rho
worst = {}
for rho, box, sizes, trials in [(2, 4, range(3, 16), 300), (3, 3, range(4, 20), 200)]:
    viol = 0; tight = 0; mx = 0
    for t in range(trials):
        sz = random.choice(list(sizes))
        pts = list(itertools.product(range(box), repeat=rho)); S = random.sample(pts, sz)
        r_ = reg(S); o = omega_cross(S); mw = MW(S)
        assert r_ <= mw + 1e-9, (S, r_, mw)        # Theorem B (sanity)
        if r_ > o + 1e-9: viol += 1
        if abs(r_ - o) < 1e-9: tight += 1
        mx = max(mx, r_/o)
    print(f"rho={rho}: {trials} random sets in [0,{box-1}]^{rho}: reg > Omega' in {viol} cases; equality in {tight}; max reg/Omega' = {mx:.3f}", flush=True)
# structured sets: boxes, simplices, cross-polytopes, diamonds
for rho in [2, 3]:
    for R in [1, 2, 3]:
        for name, S in [("box", list(itertools.product(range(R+1), repeat=rho))),
                        ("simplex", [p for p in itertools.product(range(R+1), repeat=rho) if sum(p) <= R]),
                        ("l1-ball", [p for p in itertools.product(range(-R, R+1), repeat=rho) if sum(map(abs, p)) <= R])]:
            print(f"  rho={rho} {name} R={R}: |S|={len(S)} reg={reg(S)} Omega'={omega_cross(S):.3f} MW={MW(S):.3f}", flush=True)
