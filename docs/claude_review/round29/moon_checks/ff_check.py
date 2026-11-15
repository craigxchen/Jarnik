"""Exact check of the function-field cut bound (U) and of forced collisions at degree-1 primes, over F_3.
Reuses ff_scan.py's field/polynomial code (executed with a small degree bound), then for every squarefree split N
with deg N <= 12 and every ball-cluster (fibre) of size >= 2 at c = 0, verifies:
  (i)  deg(z_x - z_y) = deg gcd + deg c_xy  with c_xy = (z_x - z_y)/gcd  (division exact);
  (ii) sum_pairs deg c_xy <= binom(M,2)*c + M*D/4 - sum_pi deg(pi) * (M - 2 S_pi)^2 / 4   [squarefree: one level];
  (iii) for l = t - a (a in F_3, inert, 4 residue classes of norm N(a)), l not dividing N: #pairs with l | c_xy
        >= pigeonhole minimum.
"""
import sys, itertools
sys.argv = ['ff_scan.py', '3', '2', '3']
src = open('ff_scan.py').read().split('best = defaultdict(int)')[0]
exec(src)
def deg(f):
    d = len(f) - 1
    while d >= 0 and f[d] == 0: d -= 1
    return d
def psub(f, g):
    n = max(len(f), len(g)); f = f + [0]*(n-len(f)); g = g + [0]*(n-len(g))
    return [ADD[a][NEG[b]] for a, b in zip(f, g)]
INV = {x: y for x in range(1, Q2) for y in range(1, Q2) if MUL[x][y] == 1}
def pdivmod(f, g):
    f = f[:deg(f)+1]; g = g[:deg(g)+1]; qd = deg(f) - deg(g)
    if qd < 0: return [0], f
    qq = [0]*(qd+1); r = f[:]; il = INV[g[-1]]
    for i in range(qd, -1, -1):
        co = MUL[r[i+len(g)-1]][il]; qq[i] = co
        for j, b in enumerate(g): r[i+j] = ADD[r[i+j]][NEG[MUL[co][b]]]
    return qq, r
def peval(f, a):
    v = 0
    for co in reversed(f): v = ADD[MUL[v][a]][co]
    return v
ok = 0; bad = 0; nclusters = 0
def rec(start, degleft, chosen):
    global ok, bad, nclusters
    if chosen:
        pis = [split[i][1] for i in chosen]; D = sum(split[i][0] for i in chosen)
        pts = [([1], ())]
        for pi in pis:
            pb = pconj(pi)
            pts = [(pmul(z, pi), e + (1,)) for z, e in pts] + [(pmul(z, pb), e + (-1,)) for z, e in pts]
        lo = int(D/2) + 1
        fib = defaultdict(list)
        for z, e in pts: fib[tuple(z[lo:D])].append((z, e))
        Nf = pmul(pts[0][0], pconj(pts[0][0]))
        for cl in fib.values():
            M = len(cl)
            if M < 2: continue
            nclusters += 1; tot = 0; coll = defaultdict(int)
            for (z1, e1), (z2, e2) in itertools.combinations(cl, 2):
                g = [1]
                for pi, a, b in zip(pis, e1, e2):
                    if a == b: g = pmul(g, pi if a == 1 else pconj(pi))
                diff = psub(z1, z2); cq, r = pdivmod(diff, g)
                assert deg(r) < 0 and deg(diff) == deg(g) + deg(cq) <= D/2
                tot += deg(cq)
                for a in range(q):
                    if peval(cq, a) == 0: coll[a] += 1
            cut = sum(split[i][0] * (M - 2*sum(1 for _, e in cl if e[j] == 1))**2 / 4 for j, i in enumerate(chosen))
            if tot <= M*D/4 - cut + 1e-9: ok += 1
            else: bad += 1
            for a in range(q):
                if peval(Nf, a) == 0: continue
                # residues z(a) in F_9^* with norm N(a): 4 classes
                qq, r = divmod(M, 4); forced = r*(qq+1)*qq//2 + (4-r)*qq*(qq-1)//2
                assert coll[a] >= forced, (coll[a], forced)
    for i in range(start, len(split)):
        d = split[i][0]
        if 2*d <= degleft:
            chosen.append(i); rec(i+1, degleft-2*d, chosen); chosen.pop()
rec(0, 12, [])
print(f"clusters checked: {nclusters}; (U) holds: {ok}; violations: {bad}; (i) and (iii) asserted for all")
