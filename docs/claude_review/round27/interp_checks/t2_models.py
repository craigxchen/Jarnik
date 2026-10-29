import sys, itertools, time
sys.path.insert(0, '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round27/interp_checks')
from shiftlib import *
from models import *
name = sys.argv[1]; Lh = int(sys.argv[2]); which = sys.argv[3] if len(sys.argv) > 3 else 'small'
if name == 'sdc7': M, pats = sdc_patterns(7)
elif name == 'syl8': M, pats = syl_patterns(3)
A, inds = box(M, pats, Lh)
d = M - 1
print(name, 'L', Lh, '#patterns', len(pats), '|A|', len(A), flush=True)
H = []
for e in range(0, 14):
    H.append(HF(A, e, d))
    if H[-1] == len(A): break
Hd = [H[0]] + [H[i] - H[i-1] for i in range(1, len(H))]
print('H_A (graded):', Hd, flush=True)
rho = max(i for i, h in enumerate(Hd) if h > 0)
if which == "pair":
    L = []
    for x in range(d):
        for s in (1, -1):
            l = [0]*d; l[x] = s; L.append(tuple(l))
    for x in range(d):
        for y in range(x+1, d):
            for s in (1, -1):
                l = [0]*d; l[x] = s; l[y] = -s; L.append(tuple(l))
elif which == "small":
    L = []
    for x in range(d):
        for s in (1, -1):
            l = [0]*d; l[x] = s; L.append(tuple(l))
    for x in range(d):
        for y in range(x+1, d):
            for s in (1, -1):
                for t in (1, -1):
                    l = [0]*d; l[x] = s; l[y] = t; L.append(tuple(l))
else:
    L = [l for l in itertools.product([-1,0,1], repeat=d) if any(l)]
print('#shifts', len(L))
for e in range(0, rho + 1):
    t0 = time.time()
    rA, racc = shift_test(A, e, d, L)
    print('e', e, 'dim M_e(A)', rA, 'dim sum M_e(A_l)', racc, 'SHIFT' if rA == racc else 'FAIL', '%.1fs' % (time.time()-t0), flush=True)
