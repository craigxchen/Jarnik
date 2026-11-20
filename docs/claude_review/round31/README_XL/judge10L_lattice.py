"""[J] Checks of the meaning referee's item (6) (column orthogonality / lattice M_S).
Certificate lattice for a cut pattern S (B rows, s primes): Lam = {a in Z^B : sum a = 0 and, for every p,
sum_{i : S_ip = 1} a_i = 0 mod p} = {a : sum a = 0, w.a = 0 mod q}, w_i = CRT_p(S_ip).
A vector a in Lam with ||a||_1 * L/2 < q is a *certificate*: a.y = 0 for every y in M_S of spread <= L
(in particular for every genuine configuration with pattern S, any residues).
We LLL-reduce Lam and count independent certificates among the reduced vectors (a lower bound for
the true number).  B - 1 independent certificates force y constant.
Patterns: (A) random S with B < s/2 rows (the regime the referee calls 'empty'); (B) the referee's
shadow patterns; (C) the five-mark genuine example (expect exactly B - 2)."""
import random
from math import isqrt, prod
from sympy import ZZ
from sympy.polys.matrices import DomainMatrix

random.seed(7)


def primes_in(lo, hi):
    return [p for p in range(lo, hi) if all(p % r for r in range(2, isqrt(p) + 1))]


def crt(res, ps):
    q = prod(ps)
    return sum(r * (q // p) * pow(q // p, -1, p) for r, p in zip(res, ps)) % q


def certificates(S, ps, C):
    B = len(S)
    q = prod(ps)
    L = isqrt(C * C * q)
    w = [crt([row[k] for k in range(len(ps))], ps) for row in S]
    u = [(w[i] - w[B - 1]) % q for i in range(B - 1)]
    K = q * 10 ** 6
    rows = []
    for i in range(B - 1):
        r = [0] * (B - 1) + [K * u[i]]
        r[i] = 1
        rows.append(r)
    rows.append([0] * (B - 1) + [K * q])
    M = DomainMatrix([[ZZ(x) for x in r] for r in rows], (B, B), ZZ).lll().to_Matrix()
    vecs = []
    for i in range(B):
        r = [int(x) for x in M.row(i)]
        if r[-1] != 0:
            continue
        c = r[:-1]
        a = c + [-sum(c)]
        vecs.append(a)
    good = [a for a in vecs if sum(abs(x) for x in a) * L < 2 * q]
    # rank of the good certificates
    if good:
        rk = DomainMatrix([[ZZ(x) for x in a] for a in good], (len(good), B), ZZ).rank()
    else:
        rk = 0
    return len(vecs), rk


out = []
pool = primes_in(50, 200)
# (A) random patterns, B < s/2
for s, B in ((14, 5), (14, 6), (16, 7), (20, 9)):
    res = []
    for t in range(6):
        ps = sorted(random.sample(pool, s))
        S = [[random.randint(0, 1) for _ in ps] for _ in range(B)]
        nv, rk = certificates(S, ps, 1)
        res.append(rk)
    out.append("random S, s=%d, B=%d (B < s/2): ranks of certificate sets over 6 trials %s (B-1=%d)" % (s, B, res, B - 1))
# (B) shadow patterns: D primes with arbitrary bits, other primes bit = [z_i = d]
for s, B in ((14, 6), (16, 8)):
    res = []
    for t in range(6):
        ps = sorted(random.sample(pool, s))
        nD = s // 2 - 1
        Z = [random.randint(0, 1) for _ in range(B)]
        Z[0], Z[1] = 0, 1
        S = [[random.randint(0, 1) for _ in range(nD)] + [Z[i]] * (s - nD) for i in range(B)]
        nv, rk = certificates(S, ps, 1)
        res.append(rk)
    out.append("shadow S, s=%d, B=%d: ranks %s (<= B-2=%d expected)" % (s, B, res, B - 2))
# (C) five-mark genuine example at C = 6
ps = [59, 67, 71, 83]
marks = [0, 11786, 11859, 23645, 28542]
S = [[0 if y % p == 0 else 1 for p in ps] for y in marks]
nv, rk = certificates(S, ps, 6)
out.append("five-mark genuine example, C=6: certificate rank %d (B-2 = 3 expected if k = 1)" % rk)
print("\n".join(out))
