# Referee's independent rerun of Evidence 2 (pair demands of the coherent G=I design with root repair
# on the sharp section-4 Paley profile, b=6). Evidence only. Own code.
import random, math, sys
sys.path.insert(0, '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round7/referee_p3Dconstruct_checks')
from paley_snf import paley
from tail_ref_lib import is_prime, leg, mul, gpow, ratio
random.seed(int(sys.argv[2]) if len(sys.argv) > 2 else 3)
q0 = int(sys.argv[1]); M = q0 + 1; b = 6; X = 3*M
H = paley(q0)
col = {a: [(1 + H[x][a+1])//2 for x in range(M)] for a in range(q0)}   # label a in F_q0 -> column; row 0 = infinity
cols = []
astar = 0
for a in range(q0):
    for i in range(b):
        c = col[a][:]
        if i == 0:                      # flip column f(x) for finite x = a (row index a+1)
            c[a+1] = 1 - c[a+1]
        if a == astar and i == 2:       # flip column f(infinity)
            c[0] = 1 - c[0]
        cols.append(c)
r = len(cols)
rows = [[cols[j][x] for j in range(r)] for x in range(M)]
assert len(set(tuple(rw) for rw in rows)) == M
for j in range(r): assert 0 < sum(cols[j]) < M
# twins check: a_{.,f(x)} - a_{.,s(x)} = +-e_x
qs = [q for q in range(3, X+1, 2) if is_prime(q)]
aux = []; n0 = X + 1
while len(aux) < r:
    if n0 % 4 == 1 and is_prime(n0):
        a = 1
        while True:
            b2 = n0 - a*a; bb = math.isqrt(b2)
            if bb*bb == b2: aux.append((a, bb)); break
            a += 1
    n0 += 1
ps = []
while len(ps) < r:
    p = random.randrange(M**4//2, M**4)
    if p % 4 == 1 and is_prime(p): ps.append(p)
npairs = M*(M-1)//2
D1 = [0.0]*npairs; DG = [0.0]*npairs
for q in qs:
    A = 1
    while q**(A+1) <= X: A += 1
    Q = q**A; chi = 1 if q % 4 == 1 else -1
    m = (q - chi)*q**(A-1); k = (m & -m).bit_length()-1; o = m >> k
    while True:
        z = (random.randrange(Q), random.randrange(Q))
        if (z[0]**2 + z[1]**2) % q == 0: continue
        f = gpow(ratio(z, Q), o, Q)
        if gpow(f, 2**(k-1), Q) != (1, 0): break
    t = []
    for j in range(r):
        Nz = aux[j][0]**2 + aux[j][1]**2
        tj = ratio(aux[j], Q)
        if leg(Nz, q) != leg(ps[j], q): tj = mul(tj, f, Q)
        t.append(tj)
    T = []
    for x in range(M):
        cur = (1, 0)
        for j in range(r):
            if rows[x][j]: cur = mul(cur, t[j], Q)
        T.append(cur)
    def inv(u):
        n = (u[0]*u[0] + u[1]*u[1]) % Q   # unit of norm 1 in T: inverse = conj
        return (u[0], (-u[1]) % Q)
    idx = 0
    lq = math.log(q)
    for x in range(M):
        Tix = inv(T[x])
        for y in range(x+1, M):
            u = mul(T[y], Tix, Q)
            u4 = gpow(u, 4, Q)
            a1 = 0; ag = 0
            for aa in range(A, 0, -1):
                Qa = q**aa
                if ag == 0 and u4[0] % Qa == 1 % Qa and u4[1] % Qa == 0: ag = aa
                if a1 == 0 and u[0] % Qa == 1 % Qa and u[1] % Qa == 0: a1 = aa
            D1[idx] += a1*lq; DG[idx] += ag*lq
            idx += 1
L = math.log(M)
print(f"q0={q0} M={M} r={r} X={X}: D1 max/mean = {max(D1)/L:.2f}/{sum(D1)/npairs/L:.3f} (units of log M); some-unit max = {max(DG)/L:.2f}")
