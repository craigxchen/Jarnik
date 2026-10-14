# gadget zonotope lower bound for two-cut blocks; compare with rho from twocut_rho_out.txt
import itertools, re, functools, sys
def gadgets():
    regs = 'PQR'
    inS1 = {'P':1,'Q':1,'R':0}; inS2 = {'P':0,'Q':1,'R':1}
    G = []
    for r in regs: G.append(((r,), 2, (inS1[r], inS2[r])))
    for r1, r2 in itertools.combinations_with_replacement(regs, 2):
        u1 = 1 if inS1[r1] != inS1[r2] else 0; u2 = 1 if inS2[r1] != inS2[r2] else 0
        G.append(((r1, r2), 2, (u1, u2)))
        h1 = 1 if (inS1[r1] or inS1[r2]) else 0; h2 = 1 if (inS2[r1] or inS2[r2]) else 0
        G.append(((r1, r2), 3, (h1, h2)))  # hexpair
    for t in itertools.combinations_with_replacement(regs, 3):
        c1 = sum(inS1[r] for r in t); c2 = sum(inS2[r] for r in t)
        G.append((t, 3, (1 if c1 in (1, 2) else 0, 1 if c2 in (1, 2) else 0)))
    return G
G = gadgets()
@functools.lru_cache(maxsize=None)
def best(a, b, c, l1, l2, i=0):
    if i == len(G): return 0
    regs, deg, (u1, u2) = G[i]
    need = {'P': regs.count('P'), 'Q': regs.count('Q'), 'R': regs.count('R')}
    res = best(a, b, c, l1, l2, i+1)
    k = 1
    while True:
        a2, b2, c2 = a - k*need['P'], b - k*need['Q'], c - k*need['R']
        m1, m2 = l1 - k*u1, l2 - k*u2
        if min(a2, b2, c2, m1, m2) < 0: break
        res = max(res, k*deg + best(a2, b2, c2, m1, m2, i+1))
        k += 1
    return res
bad = 0
for line in open('twocut_rho_out.txt'):
    m = re.match(r'a=(\d+) b=(\d+) c=(\d+) L1=(\d+) L2=(\d+) \|A\|=(\d+) rho=(\d+)', line)
    if not m: continue
    a, b, c, L1, L2, N, rho = map(int, m.groups())
    z = best(a, b, c, L1, L2)
    flag = '' if z == rho else '  <-- gap'
    if z != rho: bad += 1
    print('a=%d b=%d c=%d L1=%d L2=%d rho=%d gadget=%d%s' % (a, b, c, L1, L2, rho, z, flag))
print('cases with gadget < rho:', bad)
