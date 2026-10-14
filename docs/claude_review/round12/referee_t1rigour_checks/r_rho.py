# Referee: exact regularity rho(Lambda) computations.
# (1) near-twin sign columns: rho_A(F) for F = (1_S1, 1_S2) vs the claimed lower bound 2(js+e) used in Prop 5.5.
# (2) projection can raise rho: Lambda = {(0,0),(1,1),(2,0)}.
# (3) Lemma 6.9 values for small s (independent).
import itertools, flint, sys
from math import comb

def monos(s, e):
    out = []
    def rec(i, rem, cur):
        if i == s - 1:
            out.append(tuple(cur + [rem])); return
        for a in range(rem, -1, -1):
            rec(i + 1, rem - a, cur + [a])
    for deg in range(e + 1):
        if s == 0: out.append(()); break
        rec(0, deg, [])
    return out

def full_rank(pts, e):
    s = len(pts[0])
    ms = monos(s, e)
    if len(ms) < len(pts): return False
    M = flint.fmpz_mat([[1 if False else int(eval_mono(p, m)) for m in ms] for p in pts])
    return M.rank() == len(pts)

def eval_mono(p, m):
    v = 1
    for a, b in zip(p, m): v *= a ** b
    return v

def rho(pts):
    pts = [tuple(p) for p in set(map(tuple, pts))]
    e = 0
    while not full_rank(pts, e): e += 1
    return e

print("(2) rho({(0,0),(1,1),(2,0)}) =", rho([(0,0),(1,1),(2,0)]), " rho of its x-projection {0,1,2} =", rho([(0,),(1,),(2,)]))

print("(1) near-twin columns: Lambda = {(t1,t2): |t1|,|t2| <= j, |t1-t2| <= 2*delta} (exact image when |S1 cap S2| >= 2j+2delta)")
for j in range(1, 7):
    for delta in range(1, j + 1):
        pts = [(a, b) for a in range(-j, j+1) for b in range(-j, j+1) if abs(a-b) <= 2*delta]
        r = rho(pts)
        print(f"  j={j} delta={delta}: rho={r}  Newton bound 2j+2delta={2*j+2*delta}  width cost 4j={4*j}, alt family 2j+4delta={2*j+4*delta}"
              f"  claimed lower bounds 2(js+e): F=(2,0)->{4*j}, F'=(1,2delta)->{2*j+4*delta}")

# verify that the image is exact for a concrete box: d coordinates, S1 = I u {x}, S2 = I u {y}, |I| = 2j+2
def image(j, delta):
    I = 2*j + 2
    # coordinates: I common, delta in S1 only, delta in S2 only; others irrelevant (free cube)
    ts = set()
    for kI in range(-I, I+1):
        for k1 in range(-delta, delta+1):
            for k2 in range(-delta, delta+1):
                t1, t2 = kI + k1, kI + k2
                if abs(t1) <= j and abs(t2) <= j: ts.add((t1, t2))
    return ts
for j, delta in ((2,1),(3,1),(4,2)):
    pts = [(a, b) for a in range(-j, j+1) for b in range(-j, j+1) if abs(a-b) <= 2*delta]
    print(f"  image check j={j} delta={delta}: image == band: {image(j,delta) == set(pts)}")

print("(3) Lemma 6.9")
for s in range(2, 7):
    for t in range(0, s+1):
        L = [k for k in itertools.product((-1,0,1), repeat=s) if sum(k) == t]
        assert rho(L) == s - t, (s, t)
    for j in range(0, s+1):
        L = [k for k in itertools.product((-1,0,1), repeat=s) if abs(sum(k)) <= j]
        r = rho(L)
        assert r == s + j, (s, j, r)
print("  slices rho = s-|t| and rho(A2(s,j)) = s+j for 2<=s<=6: confirmed")
