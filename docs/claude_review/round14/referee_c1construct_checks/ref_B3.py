# Referee check of Thm B3 (independent code). L4K Thm 2 vectors a_j(q) = 2(Lam + (K+1) sigma_q(j)) + nu_j(q),
# beta_j = j, B = r, Lam = K(K+1)B + K.  F = {e_i+e_j-e_k-e_l}, unordered disjoint pairs, v ~ -v.
import itertools, random
from math import log, comb
from sympy import primerange, legendre_symbol
from fractions import Fraction

def Zr(r):
    N = {}
    for x, y in itertools.combinations(range(1, r+1), 2):
        N[x+y] = N.get(x+y, 0) + 1
    return Fraction(sum(n*n for n in N.values()), 6) - Fraction(comb(r, 2), 2)

def Fset(r):
    pairs = list(itertools.combinations(range(r), 2))
    out = []
    for P, Q in itertools.combinations(pairs, 2):
        if set(P) & set(Q): continue
        out.append((P, Q))
    return out

for r in [8, 10, 16, 24]:
    F = Fset(r); assert len(F) == r*(r-1)*(r-2)*(r-3)//8
    z = Zr(r)
    # minimise #zeros over random sigma and nu, also adversarial nu = random; check lower bound
    random.seed(r)
    K = 3; Lam = K*(K+1)*r + K
    mn = None
    for trial in range(400):
        sig = list(range(1, r+1)); random.shuffle(sig)
        nu = [random.randint(0, 1) for _ in range(r)]
        a = [2*(Lam + (K+1)*sig[j]) + nu[j] for j in range(r)]
        cnt = sum(1 for (P, Q) in F if a[P[0]] + a[P[1]] - a[Q[0]] - a[Q[1]] == 0)
        mn = cnt if mn is None else min(mn, cnt)
        assert cnt >= z
    print(f"r={r}: |F|={len(F)}, Z_r={float(z):.3f} (r^3/36={r**3/36:.2f}), min #zeros over 400 random (sigma,nu) = {mn}; Z_r/|F| * 9r/2 = {float(z/len(F))*9*r/2:.3f}")

# One full instance: K=3, r=10, primes random split <= 6000, q_T = 59, X = 60000: demand vs budget for every v in F
r, K, X, qT = 10, 3, 60000, 59
random.seed(7)
ps = sorted(random.sample([p for p in primerange(1000, 6000) if p % 4 == 1], r))
Lam = K*(K+1)*r + K
F = Fset(r)
dem = {f: 0.0 for f in F}
for q in primerange(qT+1, X+1):
    if q in ps: continue
    Aq = 1
    while q**(Aq+1) <= X: Aq += 1
    sig = list(range(1, r+1)); random.shuffle(sig)
    nu = [1 if legendre_symbol(p, q) == -1 else 0 for p in ps]
    a = [2*(Lam + (K+1)*sig[j]) + nu[j] for j in range(r)]
    for (P, Q) in F:
        if a[P[0]] + a[P[1]] - a[Q[0]] - a[Q[1]] == 0:
            dem[(P, Q)] += Aq*log(q)
bud = {(P, Q): 0.5*sum(log(ps[i]) for i in P+Q) for (P, Q) in F}
over = sum(1 for f in F if dem[f] > bud[f])
worst = min(F, key=lambda f: dem[f]-bud[f])
print(f"instance K=3 r=10 primes {ps}, q in ({qT},{X}]: monomials of F over unit-0 budget: {over}/{len(F)}; "
      f"min demand {min(dem.values()):.1f}, max budget {max(bud.values()):.2f}")
