"""Translate an integral-cofactor configuration z_i (= (A_i - A_1)/G) by a Gaussian
integer c so that every z_i + c has EXACT exponent row: divisible by f_ip exactly once,
not divisible by conj f_ip, for every p | q. CRT over the 2s Gaussian primes pi_p, conj pi_p
(and pi_p^2 to keep exponent exactly 1). Usage: python3 rows.py P i j k ..."""
import math, pickle, sys, random
from marks import gmul

P = int(sys.argv[1]); idx = [int(a) for a in sys.argv[2:]]
ps, pis, q, marks = pickle.load(open(f"marks_{P}.pkl", "rb"))
s = len(ps)
G = 0
for a in idx:
    for b in idx:
        if a < b:
            dx = marks[b][0] - marks[a][0]; dy = marks[b][1] - marks[a][1]
            G = math.gcd(G, math.gcd(dx, dy) // 2)
z = {r: ((marks[r][0] - marks[idx[0]][0]) // G, (marks[r][1] - marks[idx[0]][1]) // G) for r in idx}

def conj(w): return (w[0], -w[1])
def val(w, pi):
    """pi-adic valuation of Gaussian integer w (w != 0)."""
    n = pi[0] ** 2 + pi[1] ** 2; v = 0
    while True:
        num = gmul(w, conj(pi))
        if num[0] % n or num[1] % n:
            return v
        w = (num[0] // n, num[1] // n); v += 1

# Work in Z[i]/p^2 ~ via the rational representation: Z[i]/(pi^2) = Z/p^2 with i -> root.
def root_mod(pi, m):
    # i = -a/b mod pi^k: pi = a + b i, so i = -a * b^{-1} mod (pi^k) in Z/p^k
    a, b = pi
    p = a * a + b * b
    r = (-a * pow(b, -1, p)) % p          # i = r mod pi
    # Hensel lift of r^2 + 1 = 0 from p to p^2 (m = p^2)
    r = (r - (r * r + 1) * pow(2 * r, -1, m)) % m
    assert (r * r + 1) % m == 0 and (a + b * r) % p == 0
    return r
residues = []  # list of (modulus p^2, root r, target residue) for c = x + y i -> x + y r
for k, (p, pi) in enumerate(zip(ps, pis)):
    for g in (pi, conj(pi)):
        m = p * p; r = root_mod(g, m)
        members = [t for t in idx if (((t >> k) & 1) == 1) == (g == pi)]
        if not members:
            # no mark uses this conjugate: choose c mod p so that no z_t + c is divisible by g
            badp = {(-(z[t][0] + z[t][1] * r)) % p for t in idx}
            cand = next(h for h in range(p) if h not in badp)
            residues.append((m, r, cand)); continue
        base = [(-(z[t][0] + z[t][1] * r)) % p for t in members]
        assert len(set(base)) == 1, "consistency mod f_ip fails"
        # choose lift mod p^2 avoiding z_t + c = 0 mod g^2
        bad = {(-(z[t][0] + z[t][1] * r)) % m for t in members}
        for h in range(p):
            cand = (base[0] + p * h) % m
            if cand not in bad:
                residues.append((m, r, cand)); break
# CRT: x + y r = cand mod p^2 for each (p^2, r). Take y = 0 -> x = cand mod each modulus? moduli
# for pi and conj pi are both p^2 with different roots; solve with y as well: two congruences mod p^2
# x + y r1 = c1, x + y r2 = c2  ->  y = (c1 - c2)/(r1 - r2), r1 - r2 invertible mod p (r1 = -r2).
X = 0; Y = 0; MOD = 1
byp = {}
for (m, r, c) in residues:
    byp.setdefault(m, []).append((r, c))
for m, lst in byp.items():
    if len(lst) == 2:
        (r1, c1), (r2, c2) = lst
        y = ((c1 - c2) * pow(r1 - r2, -1, m)) % m; x = (c1 - y * r1) % m
    else:
        (r1, c1), = lst; y = 0; x = c1
    # CRT combine
    X = X + MOD * (((x - X) * pow(MOD, -1, m)) % m)
    Y = Y + MOD * (((y - Y) * pow(MOD, -1, m)) % m)
    MOD *= m
c = (X, Y)
ok = True
for t in idx:
    w = (z[t][0] + c[0], z[t][1] + c[1])
    row = []
    for k, pi in enumerate(pis):
        want = pi if (t >> k) & 1 else conj(pi)
        other = conj(pi) if (t >> k) & 1 else pi
        e1, e2 = val(w, want), val(w, other)
        ok &= (e1 == 1 and e2 == 0)
        row.append((e1, e2))
    print(t, "z+c =", w, "norm/q =", (w[0] ** 2 + w[1] ** 2) / q, "rows exact:", all(e == (1, 0) for e in row))
print("translate c =", c, "| all exponent rows exact:", ok)
