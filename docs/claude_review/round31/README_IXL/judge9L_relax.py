# Judge spot check (d): IX.4 relaxation at m = 11 and m = 19. Own code, exact integers.
# Paley cuts S = K - I, K_iv = chi(v - i); two primes per column; edges d_ij = h_ij (w_j - w_i).
# Checks: S^T S = (m+1)I - J; Sidon pair sums; gcd(d_ij, q) = h_ij; triangle failures;
# pair-level signed residue test: a genuine two-class configuration has d_ij^2 mod p constant over separated pairs.
import math, itertools
out = []
def say(*a):
    s = " ".join(str(x) for x in a); print(s); out.append(s)
def chi(a, m):
    a %= m
    if a == 0: return 0
    return 1 if pow(a, (m-1)//2, m) == 1 else -1
def primes(lo, hi):
    return [p for p in range(lo, hi) if all(p % d for d in range(2, int(p**.5)+1))]
for m in (11, 19):
    S = [[(chi(v-i, m) if v != i else -1) for v in range(m)] for i in range(m)]
    StS = [[sum(S[k][i]*S[k][j] for k in range(m)) for j in range(m)] for i in range(m)]
    okgram = all(StS[i][j] == ((m+1) if i == j else 0) - 1 for i in range(m) for j in range(m))
    P = 4000
    ps = [p for p in primes(P, 2*P) if p % 4 == 1][:2*m]
    cols = [(ps[2*v], ps[2*v+1]) for v in range(m)]
    q = math.prod(ps)
    bits = [[(1+S[i][v])//2 for v in range(m)] for i in range(m)]
    sig = [tuple(b for b in bits[i] for _ in (0, 1)) for i in range(m)]
    sums = {tuple(a ^ b for a, b in zip(sig[i], sig[j])) for i, j in itertools.combinations(range(m), 2)}
    h = {(i, j): math.prod(cols[v][0]*cols[v][1] for v in range(m) if bits[i][v] == bits[j][v])
         for i in range(m) for j in range(m) if i != j}
    k = P // (4*m)
    w = [k*i + (1 if i == m-1 else 0) for i in range(m)]
    d = {(i, j): h[(i, j)]*(w[j]-w[i]) for (i, j) in h}
    okgcd = all(math.gcd(d[(i, j)], q) == h[(i, j)] for (i, j) in d)
    fails = sum(1 for i, j, l in itertools.combinations(range(m), 3) if d[(i, j)] + d[(j, l)] != d[(i, l)])
    distinct = []
    for v in range(m):
        for p in cols[v]:
            sq = {d[(i, j)]**2 % p for i, j in itertools.combinations(range(m), 2) if bits[i][v] != bits[j][v]}
            npairs = sum(1 for i, j in itertools.combinations(range(m), 2) if bits[i][v] != bits[j][v])
            distinct.append((len(sq), npairs))
    say(f"m={m}: Gram identity {okgram}; distinct pair sums {len(sums)} of {m*(m-1)//2}; gcd(d,q)=h {okgcd}; "
        f"failed triangles {fails} of {math.comb(m,3)}; max|d|/sqrt q = {max(abs(x) for x in d.values())/math.isqrt(q):.3f}")
    say(f"  d_ij^2 mod p over separated pairs: (#distinct values, #pairs) range "
        f"{min(distinct)} .. {max(distinct)}; a genuine configuration has 1 value at every p")
open(__file__.replace('.py', '.out'), 'w').write("\n".join(out) + "\n")
