"""Referee exploration: do the PAIR inequalities alone bound the flip fraction W_F/W in fully flipped
Hadamard-core profiles?  (Relevant to the 'heavy-flip regime W_F >= W/2' of flipped.md Sec. 7.)

Model (one flip per row at distinct columns, assignment x -> a(x), capacity <= B):
  variables  U_a >= 0  (unflipped weight of label a),   l_x >= 0 (log of the flipped prime of row x);
  W_a = U_a + sum_{a(x)=a} l_x,  W = sum_a U_a + sum_x l_x,  W_F = sum_x l_x.
  Pair inequality (walsh-extraction Sec. 1, G_xy <= 4 log C - D), x != y:
     G_xy = sum_a W_a H(x,a)H(y,a) + 2 eps_xy l_x + 2 eps_yx l_y,   eps_xy = -H(x,a(x)) H(y,a(x)).
  Asymptotic (scale-free) relaxation: G_xy <= 0.  (Exact constraint is G_xy <= 4 log C - D, which is
  o(W) at fixed C.)
LP: maximise W_F subject to G_xy <= 0 (all pairs), W <= 1, variables >= 0.  The optimum is the largest
flip fraction W_F/W compatible with all pair inequalities (in the limit W -> infinity at fixed C).
Exact rational simplex (Bland's rule).  This is a RELAXATION: it ignores integrality/primality of weights.
"""
from fractions import Fraction as Fr
import random, sys

def legendre(a, q):
    a %= q
    if a == 0: return 0
    return 1 if pow(a, (q-1)//2, q) == 1 else -1
def paley(q):
    M = q + 1
    return [[1]*M] + [[1] + [-legendre(a - x, q) - (a == x) for a in range(q)] for x in range(q)]
def sylv(t):
    M = 2**t
    return [[(-1)**bin(x & a).count('1') for a in range(M)] for x in range(M)]

def simplex_max(A, b, c):
    """max c.z s.t. A z <= b, z >= 0, with b >= 0 (origin feasible). Exact Fractions, Bland's rule."""
    m, n = len(A), len(c)
    T = [[Fr(v) for v in A[i]] + [Fr(int(i == k)) for k in range(m)] + [Fr(b[i])] for i in range(m)]
    obj = [Fr(-v) for v in c] + [Fr(0)]*m + [Fr(0)]
    basis = [n + i for i in range(m)]
    while True:
        e = next((j for j in range(n + m) if obj[j] < 0), None)
        if e is None: break
        best, r = None, None
        for i in range(m):
            if T[i][e] > 0:
                ratio = T[i][-1] / T[i][e]
                if best is None or ratio < best or (ratio == best and basis[i] < basis[r]):
                    best, r = ratio, i
        if r is None: return None, None  # unbounded
        piv = T[r][e]
        T[r] = [v / piv for v in T[r]]
        for i in range(m):
            if i != r and T[i][e] != 0:
                f = T[i][e]; T[i] = [vi - f*vr for vi, vr in zip(T[i], T[r])]
        if obj[e] != 0:
            f = obj[e]; obj = [vo - f*vr for vo, vr in zip(obj, T[r])]
        basis[r] = e
    z = [Fr(0)]*(n + m)
    for i, bi in enumerate(basis): z[bi] = T[i][-1]
    return obj[-1], z[:n]

def lp_for(H, assign):
    M = len(H); A_ = list(range(1, M))
    nv = (M - 1) + M            # U_a (index a-1), l_x (index M-1+x)
    rows, rhs = [], []
    for x in range(M):
        for y in range(x + 1, M):
            coef = [0]*nv
            for a in A_:
                h = H[x][a]*H[y][a]
                coef[a - 1] += h
                for x2 in range(M):
                    if assign[x2] == a: coef[M - 1 + x2] += h
            exy = -H[x][assign[x]]*H[y][assign[x]]
            eyx = -H[y][assign[y]]*H[x][assign[y]]
            coef[M - 1 + x] += 2*exy; coef[M - 1 + y] += 2*eyx
            rows.append(coef); rhs.append(0)
    rows.append([1]*nv); rhs.append(1)
    c = [0]*(M - 1) + [1]*M
    val, z = simplex_max(rows, rhs, c)
    return val, z

def random_assign(M, B, rng):
    while True:
        lab = []
        cnt = {}
        ok = True
        for x in range(M):
            choices = [a for a in range(1, M) if cnt.get(a, 0) < B]
            a = rng.choice(choices); lab.append(a); cnt[a] = cnt.get(a, 0) + 1
        return lab

if __name__ == '__main__':
    rng = random.Random(5)
    cases = [('Walsh-8', sylv(3)), ('Paley-12', paley(11)), ('Walsh-16', sylv(4)), ('Paley-20', paley(19))]
    for name, H in cases:
        M = len(H)
        res = []
        for trial in range(6 if M <= 12 else 3):
            if trial == 0:
                assign = [1] + list(range(1, M))          # canonical capacity-two
            else:
                labels = list(range(1, M)) + [rng.randrange(1, M)]; rng.shuffle(labels); assign = labels[:M]
            val, z = lp_for(H, assign)
            U = sum(z[:M - 1]); L = sum(z[M - 1:])
            res.append(val)
            print(f"{name} assign#{trial}: max W_F/W under all pair inequalities = {val} = {float(val):.4f}"
                  f"   (min label U_a/W at optimum = {float(min(z[:M-1])):.4f})", flush=True)
        print(f"{name}: range of LP optima {float(min(res)):.4f} .. {float(max(res)):.4f}")
