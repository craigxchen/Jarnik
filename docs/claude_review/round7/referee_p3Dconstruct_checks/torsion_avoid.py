# M2(b): labels uniform on (parity class) minus 4-torsion: single-column hit prob 0;
# pair (same parity) hit prob exactly 4/(m-4) (k=2) / 8/(m-8) (k>=3), >= p(q,0). Exact enumeration.
from fractions import Fraction as F
def check(m):
    k = (m & -m).bit_length()-1
    T4 = {0, m//4, m//2, 3*m//4}
    out = []
    for nu in (0, 1):
        S = [x for x in range(m) if x % 2 == nu and x not in T4]
        if not S: continue
        single = F(sum(1 for x in S if x in T4), len(S))
        pair = F(sum(1 for x in S for y in S if (x - y) % m in T4), len(S)**2)
        p0 = F(4, m) if k == 2 else F(8, m)
        out.append((nu, single, pair, pair >= p0))
    return k, out
for m in [12, 20, 28, 44, 48, 64, 100, 144, 196, 400]:
    k, out = check(m)
    print(m, "k=", k, [(nu, str(s), str(p), ok) for nu, s, p, ok in out])
