# Referee's variant: V = span(1,x2,x3,x4) (dim 4), five hyperplanes H not containing 1, F_H(Y)=P_H(Y)/Y (deg 7).
# Symbolic contact check and Jacobian rank of (x2,x3,x4,Y) -> [F_H(Y)] at random points of F_{2^20}.
import itertools, random, sys
from gf import *
random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 2)
n = 4  # variables x2,x3,x4,Y -> indices 0,1,2,3
one = mp_const1(n); xv = [mp_var(i, n) for i in range(3)]; Yv = mp_var(3, n)
def elem(c):  # c = (c1, c2, c3, c4) coefficient on (1,x2,x3,x4)
    s = set()
    if c[0]: s = mp_add(s, one)
    for i in range(3):
        if c[i + 1]: s = mp_add(s, xv[i])
    return s
for name, sel in (("a", [(0,0,0),(1,0,0),(0,1,0),(0,0,1),(1,1,1)]), ("b", [(0,0,0),(1,0,0),(0,1,0),(1,1,0),(0,0,1)]),
                  ("all8", list(itertools.product([0,1], repeat=3)))):
    Fs = []
    for lv in sel:
        r = mp_const1(n)
        for c in itertools.product([0, 1], repeat=4):
            val = c[0] ^ (c[1] & lv[0]) ^ (c[2] & lv[1]) ^ (c[3] & lv[2])
            if val == 0 and any(c): r = mp_mul(r, mp_add(Yv, elem(c)))
        Fs.append(r)
    ydeg = {e[3] for F_ in Fs for e in F_}
    mx = max(max(e[3] for e in mp_add(Fs[i], Fs[j])) for i in range(len(Fs)) for j in range(i+1, len(Fs)))
    F = GF2k(20)
    dF = [[mp_diff(Fs[j], k) for k in range(n)] for j in range(len(Fs))]
    ranks = []
    for trial in range(4):
        pt = [F.rand() for _ in range(n)]
        M = [[mp_eval(F, dF[j][k], pt) for k in range(n)] + [mp_eval(F, Fs[j], pt)] for j in range(len(Fs))]
        ranks.append(rank_gf(F, M) - 1)
    print(f"selection {name}: Y-exponents {sorted(ydeg)}; max deg_Y(F_i-F_j) = {mx} vs 7;"
          f" projective Jacobian rank at random points: {ranks}")
