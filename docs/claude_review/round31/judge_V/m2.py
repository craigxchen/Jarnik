# Exact M_2(D,0): largest fibre of monic z in F_4[t], deg D, under z -> (norm, coefficients of degree > floor(D/2)).
import sys
from collections import Counter
from f4 import mul2
def M2(D):
    h = D // 2
    cnt = Counter()
    # monic: A has t^D coefficient 1, B has degree < D
    for A0 in range(1 << D):
        A = (1 << D) | A0
        AA = mul2(A, A)
        for B in range(1 << D):
            N = AA ^ mul2(A, B) ^ mul2(B, B)
            cnt[(N, A >> (h + 1), B >> (h + 1))] += 1
    mx = max(cnt.values())
    return mx, sum(1 for v in cnt.values() if v == mx)
for D in range(1, int(sys.argv[1]) + 1):
    print(D, *M2(D), flush=True)
