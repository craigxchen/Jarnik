"""Test C_bal:  reg(S) * pi_bal <= E_bal[w_P(S)]   (pi_bal = Pr[uniform balanced cut separates a pair])
and Conjecture B: reg(S) <= 2 E_iid[w_P(S)], on random finite S (reg exact mod 2 primes = rigorous
upper bound; the inequality is checked against that upper bound, so a PASS is rigorous)."""
import sys, random, itertools
from fractions import Fraction as Fr
from clib import reg2, proj, widths, popc
from gens import orbit
random.seed(int(sys.argv[3]))
M = int(sys.argv[1]); trials = int(sys.argv[2]); R = int(sys.argv[4])
q = M // 2
def Ebal(w):
    # uniform over cuts of size q (and q+1 for odd M, same up to complement)
    vals = [v for P, v in w.items() if min(popc(P), M - popc(P)) == q]
    return Fr(sum(vals), len(vals))
def Eiid(w):
    return Fr(2 * sum(w.values()), 2**M)
pib = Fr(2 * q * (M - q), M * (M - 1))
worstb = None; worsti = None; cnt = 0
for t in range(trials):
    mode = random.choice(["rand", "orb", "box", "lattice"])
    if mode == "rand":
        S = set()
        for _ in range(random.randint(2, 40)):
            k = [random.randint(-R, R) for _ in range(M - 1)]; S.add(tuple(k) + (-sum(k),))
    elif mode == "orb":
        S = set()
        for _ in range(random.randint(1, 2)):
            v = [random.randint(-R, R) for _ in range(M - 1)]; v.append(-sum(v)); S |= set(orbit(tuple(v)))
    elif mode == "box":
        # lattice points of a random parallelepiped spanned by small integer sum-zero vectors
        gens = []
        for _ in range(random.randint(1, M - 1)):
            v = [random.randint(-2, 2) for _ in range(M - 1)]; v.append(-sum(v)); gens.append(v)
        lens = [random.randint(1, 3) for _ in gens]
        S = set()
        for coeffs in itertools.product(*[range(L + 1) for L in lens]):
            S.add(tuple(sum(c * g[i] for c, g in zip(coeffs, gens)) for i in range(M)))
    else:
        # random cut polytope lattice points
        U = {}
        for r in range(1, M):
            for P in itertools.combinations(range(M), r):
                U[P] = random.randint(0, R)
        S = set()
        for k in itertools.product(range(-R, R + 1), repeat=M - 1):
            kk = k + (-sum(k),)
            if all(sum(kk[j] for j in P) <= U[P] for P in U): S.add(kk)
    S = sorted(S)
    if len(S) < 2 or len(S) > 1500: continue
    cnt += 1
    r = reg2(proj(S)); w = widths(S, M)
    mb = Ebal(w) - r * pib; mi = 2 * Eiid(w) - r
    if worstb is None or mb < worstb[0]: worstb = (mb, mode, r, Ebal(w), S if len(S) < 30 else len(S))
    if worsti is None or mi < worsti[0]: worsti = (mi, mode, r, 2 * Eiid(w), S if len(S) < 30 else len(S))
print(f"M={M} sets={cnt}: min (E_bal w - pi_bal reg) = {worstb[0]} [{worstb[1]}, reg={worstb[2]}, Ebal={worstb[3]}] ; "
      f"min (2E_iid w - reg) = {worsti[0]} [{worsti[1]}, reg={worsti[2]}]")
print("  worst-bal set:", worstb[4])
