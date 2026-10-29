import sys, itertools, time
sys.path.insert(0, '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round27/interp_checks')
from shiftlib import *
from models import *
name = sys.argv[1]; Lh = int(sys.argv[2]); e = int(sys.argv[3])
if name == 'sdc7': M, pats = sdc_patterns(7)
elif name == 'syl8': M, pats = syl_patterns(3)
A, inds = box(M, pats, Lh)
d = M - 1
Aset = set(A)
# all shifts in A - A
D = set()
for a in A:
    for b in A:
        D.add(tuple(x - y for x, y in zip(a, b)))
D.discard(tuple([0]*d))
# keep shifts with A_l large enough to matter, sorted by |A_l| desc
sizes = []
for l in D:
    s = sum(1 for k in A if tuple(ki + li for ki, li in zip(k, l)) in Aset)
    sizes.append((s, l))
sizes.sort(reverse=True)
print(name, 'L', Lh, '|A|', len(A), '#diff shifts', len(D), 'largest |A_l|', sizes[:5][0][0], flush=True)
# only shifts with |A_l| >= dim Sym^{<e}... any A_l with M_e(A_l)!=0 must have |A_l| > HF_{A_l}(e-1)
Ls = [l for s, l in sizes if s >= 10]
t0 = time.time()
rA, racc = shift_test(A, e, d, Ls)
print('e', e, 'dim M_e(A)', rA, 'dim sum over', len(Ls), 'shifts', racc, 'SHIFT' if rA == racc else 'FAIL', '%.1fs' % (time.time()-t0), flush=True)
