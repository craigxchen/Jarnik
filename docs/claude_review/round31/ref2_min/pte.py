# Integer vectors m in Z^8 with sum m = sum m x = sum m x^2 = 0 (roots x+i: Im rho^k gives 1, 2x, 3x^2-1),
# l1 <= 16, enumerated by meet in the middle on coordinates 0-3 | 4-7. Independent of Astra's multiset code.
import itertools
from collections import defaultdict, Counter
X = [1,5,8,10,12,16,22,23]
# check moment vectors of x+i: Im (x+i)^k
def im(z): return z.imag
h = [(1, 2*x, 3*x*x-1) for x in X]
for x,hh in zip(X,h):
    z=complex(x,1); assert all(abs(im(z**k)-hh[k-1])<1e-9 for k in (1,2,3))
W = 16
def half(idx):
    out = defaultdict(list)
    def rec(pos, m, w):
        if pos == 4:
            s = tuple(sum(mi*h[i][c] for mi,i in zip(m,idx)) for c in range(3))
            out[s].append((tuple(m), w)); return
        for v in range(-(W-w), W-w+1):
            rec(pos+1, m+[v], w+abs(v))
    rec(0, [], 0); return out
A = half([0,1,2,3]); B = half([4,5,6,7])
cnt = Counter(); vecs = defaultdict(list)
for s, la in A.items():
    t = tuple(-c for c in s)
    if t in B:
        for ma, wa in la:
            for mb, wb in B[t]:
                if wa+wb <= W and wa+wb > 0:
                    cnt[wa+wb] += 1; vecs[wa+wb].append(ma+mb)
print('number of nonzero kernel vectors by l1 weight (both signs):', sorted(cnt.items()))
for m in vecs[14]:
    if next(x for x in m if x) > 0: print(m)
