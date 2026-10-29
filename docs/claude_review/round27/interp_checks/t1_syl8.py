import sys, itertools, time
sys.path.insert(0, '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round27/interp_checks')
from shiftlib import *
# Sylvester n=8: points x in F_2^3, base point n = 0; columns u in F_2^3\0, beta_u(x)=(-1)^{u.x}
pts = [x for x in itertools.product([0,1], repeat=3)]
base = (0,0,0)
others = [x for x in pts if x != base]
d = len(others)
cols = [u for u in itertools.product([0,1], repeat=3) if any(u)]
def dot(u,x): return sum(a*b for a,b in zip(u,x)) % 2
B = [[((-1)**dot(u,x)) - 1 for x in others] for u in cols]   # b_{u,x} = beta_u(x)-beta_u(base)
j = int(sys.argv[1]) if len(sys.argv) > 1 else 1
w = 2*j
A = [k for k in itertools.product([-1,0,1], repeat=d) if all(abs(sum(b*kk for b,kk in zip(row,k))) <= w for row in B)]
print('d', d, 'j', j, '|A|', len(A), 'of', 3**d)
# shifts: unit vectors and pair differences/sums
L = []
for x in range(d):
    for s in (1,-1):
        l = [0]*d; l[x] = s; L.append(tuple(l))
for x in range(d):
    for y in range(x+1, d):
        for s in (1,-1):
            for t in (1,-1):
                l = [0]*d; l[x] = s; l[y] = t; L.append(tuple(l))
emax = int(sys.argv[2]) if len(sys.argv) > 2 else 8
for e in range(0, emax+1):
    t0 = time.time()
    rA, racc = shift_test(A, e, d, L)
    print('e', e, 'dim M_e(A)', rA, 'dim sum M_e(A_l)', racc, 'SHIFT' if rA == racc else 'FAIL', '%.1fs' % (time.time()-t0), flush=True)
    if rA == 0: break
