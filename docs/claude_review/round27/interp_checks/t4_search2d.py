"""Search small 2-D point sets S for: shift-criterion reach m_shift vs max binomial-product order zeta(S).
Evidence-level search (numeric complex arithmetic for zeta; mod-P ranks for the criterion)."""
import sys, itertools, random, cmath, math, time
sys.path.insert(0, '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round27/interp_checks')
from shiftlib import *

def graded_H(S, d):
    H = []
    e = 0
    while True:
        H.append(HF(S, e, d))
        if H[-1] == len(S): break
        e += 1
    return [H[0]] + [H[i] - H[i-1] for i in range(1, len(H))]

def shift_reach(S, d):
    Sset = set(S)
    L = set(tuple(a[i]-b[i] for i in range(d)) for a in S for b in S)
    L.discard(tuple([0]*d))
    L = sorted(L)
    Hd = graded_H(S, d)
    rho = len(Hd) - 1
    m = 0
    for e in range(0, rho + 1):
        rA, racc = shift_test(S, e, d, L)
        if rA != racc:
            return e + 1 if False else e, rho  # SHIFT_e fails: criterion certifies orders m <= e
    return rho + 1, rho

# binomial products: factors (v, zeta), v primitive with first nonzero coord > 0
def factors(N, K=6):
    fs = []
    for a in range(0, N + 1):
        for b in range(-N, N + 1):
            if (a, b) == (0, 0): continue
            if a == 0 and b < 0: continue
            if math.gcd(a, abs(b)) != 1: continue
            for k in range(K):
                fs.append(((a, b), cmath.exp(2j * math.pi * k / K), k == 0))
    return fs

def polymul(P, f):
    (v, z, _) = f
    Q = {}
    for m, c in P.items():
        m2 = (m[0] + v[0], m[1] + v[1])
        Q[m2] = Q.get(m2, 0) + c
        Q[m] = Q.get(m, 0) - z * c
    return {m: c for m, c in Q.items() if abs(c) > 1e-9}

def fits(P, Sset, bbox):
    # is there a translate a0 with supp(P)+a0 subset S?
    ms = list(P.keys())
    m0 = ms[0]
    for s in Sset:
        a0 = (s[0] - m0[0], s[1] - m0[1])
        if all((m[0] + a0[0], m[1] + a0[1]) in Sset for m in ms):
            return True
    return False

def zeta(S, maxt=12):
    Sset = set(S)
    xs = [s[0] for s in S]; ys = [s[1] for s in S]
    N = max(max(xs) - min(xs), max(ys) - min(ys))
    fs = factors(N)
    best = 0
    # DFS over nondecreasing factor index
    sys.setrecursionlimit(10000)
    def width(P):
        a = [m[0] for m in P]; b = [m[1] for m in P]
        return max(a) - min(a), max(b) - min(b)
    def dfs(P, start, order):
        nonlocal best
        if order > best: best = order
        for i in range(start, len(fs)):
            Q = polymul(P, fs[i])
            wx, wy = width(Q)
            if wx > max(xs) - min(xs) or wy > max(ys) - min(ys): continue
            if not fits(Q, Sset, None): continue
            dfs(Q, i, order + (1 if fs[i][2] else 0))
    dfs({(0, 0): 1.0}, 0, 0)
    return best

if __name__ == '__main__':
    random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
    G = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    trials = int(sys.argv[3]) if len(sys.argv) > 3 else 30
    grid = [(a, b) for a in range(G + 1) for b in range(G + 1)]
    for t in range(trials):
        k = random.randint(len(grid) // 2, len(grid) - 1)
        S = sorted(random.sample(grid, k))
        t0 = time.time()
        ms, rho = shift_reach(S, 2)
        z = zeta(S)
        flag = '  <-- criterion beyond zeta' if ms > z else ''
        print(t, '|S|', len(S), 'rho', rho, 'shift certifies m<=', ms, 'zeta', z, '%.1fs' % (time.time() - t0), flag, S if flag else '', flush=True)
