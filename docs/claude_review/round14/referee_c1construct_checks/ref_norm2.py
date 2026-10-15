# Referee evidence: on {5,13,17,29} with genuine labels at all odd q <= X (X large), do norm-2 relations
# e_j +- e_k with the other angles FIXED at their genuine values pin the norm-1-unpinned columns 13, 17?
# Free set of phi_j in [0,pi) = complement of arcs (s pi/4 - <v',theta> +- asin(nu_s m_s/|Z_v|)), 50 digits.
import mpmath as mp
from math import isqrt
mp.mp.dps = 50
def gauss(p):
    for b in range(2, isqrt(p)+1, 2):
        a2 = p-b*b; a = isqrt(a2)
        if a*a == a2 and a % 2: return a, b
ps = [5, 13, 17, 29]
G = {p: complex(*gauss(p)) for p in ps}
def nfree(n):
    n = abs(int(n))
    while n % 2 == 0: n //= 2
    for pl in ps:
        while n % pl == 0: n //= pl
    return n
def free_set(j, others):
    # others: list of (k, sign) ; include norm-1 (empty) too
    ivs = []
    rels = [None] + [(k, sg) for k in ps if k != j for sg in (1, -1)] if others else [None]
    for rel in rels:
        a, b = gauss(j); Z = (a, b); shift = mp.mpf(0)
        if rel:
            k, sg = rel; c, d = gauss(k)
            if sg == -1: d = -d
            Z = (a*c - b*d, a*d + b*c); shift = sg*mp.atan2(gauss(k)[1], gauss(k)[0])
        X_, Y_ = Z; absZ = mp.sqrt(X_*X_ + Y_*Y_)
        D = [Y_, X_-Y_, X_, X_+Y_]
        for s in range(4):
            nu = 1 if s % 2 == 0 else 1/mp.sqrt(2)
            r = mp.asin(min(1, nu*nfree(D[s])/absZ)); c0 = s*mp.pi/4 - shift
            for sh in range(-3, 4):
                L, H = max(c0-r+sh*mp.pi, mp.mpf(0)), min(c0+r+sh*mp.pi, mp.pi)
                if L < H: ivs.append((L, H))
    ivs.sort(); cur = mp.mpf(0); free = []
    for L, H in ivs:
        if L > cur + mp.mpf(10)**-40: free.append((float(cur), float(L)))
        cur = max(cur, H)
    if cur < mp.pi - mp.mpf(10)**-40: free.append((float(cur), float(mp.pi)))
    return free
for j in [13, 17]:
    print(f"p={j}: theta={float(mp.atan2(gauss(j)[1], gauss(j)[0])):.4f}; norm-1 free set {free_set(j, False)}; "
          f"with norm-2 relations to the other columns at genuine angles: {free_set(j, True)}")
