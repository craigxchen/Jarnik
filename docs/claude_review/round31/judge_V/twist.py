# lower's GAP: unit-twisted differences at T = 0 (Y = theta). Check F_i/theta - F_j/theta^2 = A*O exactly,
# and deg O = deg P_J(1) - a1 with a1 in {1,2}.
import itertools
from f4 import *
def par(x): return bin(x).count('1') & 1
TH2 = (1, 1)           # theta^2 = theta + 1 ; theta^{-1} = theta^2, theta^{-2} = theta
def run(d):
    Q = 2 ** d; V = range(Q); forms = [c for c in range(Q) if c & 1]
    lin = lambda v: (v, 1)                      # theta + v
    cnt = 0; losses = set()
    for c1, c2 in itertools.permutations(forms, 2):
        Si = [v for v in V if par(c1 & v) == 0]            # side of c1 containing 0
        Sj = [v for v in V if par(c2 & v) == 1]            # side of c2 containing 1
        Fi = ONE; Fj = ONE
        for v in Si: Fi = mul4(Fi, lin(v))
        for v in Sj: Fj = mul4(Fj, lin(v))
        tw = add4(mul4(Fi, TH2), mul4(Fj, THETA))          # F_i/theta + F_j/theta^2 (char 2)
        A = ONE
        for v in set(Si) & set(Sj): A = mul4(A, lin(v))
        J = [v for v in V if par(c1 & v) == 0 and par(c2 & v) == 0]
        # P_J(X) as list of coefficients of X^(2^k): build polynomial in X over F_2[t]
        P = [1]                                              # coefficients by X-degree
        for v in J:                                          # multiply by (X + v)
            newP = [0] * (len(P) + 1)
            for i, a in enumerate(P):
                newP[i + 1] ^= a; newP[i] ^= mul2(a, v)
            P = newP
        ks = [i for i, a in enumerate(P) if a]
        assert all((i & (i - 1)) == 0 for i in ks), "not additive"
        O = 0
        for k in range(20):
            if (1 << k) < len(P) and k % 2 == 1: O ^= P[1 << k]
        assert tw == mul4(A, (O, 0)), (d, c1, c2)
        PJ1 = 1
        for v in J: PJ1 = mul2(PJ1, 1 ^ v)
        losses.add(deg2(PJ1) - deg2(O)); cnt += 1
    print(f"d={d}: twisted identity exact on {cnt} ordered pairs; deg P_J(1) - deg O in {sorted(losses)}")
for d in (3, 4, 5, 6): run(d)
