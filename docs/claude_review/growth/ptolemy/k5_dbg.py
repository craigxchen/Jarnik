from k5_family_search import roots, good, coprime_class
import collections
cnt = collections.Counter()
ex = None
for lam in range(-6, 7):
    if lam in (0, -1, 1, -2): continue
    for a in range(-40, 41, 4):
        for b in range(-40, 41, 4):
            for c in range(-40, 41, 4):
                R = roots(a, b, c, lam)
                if R is None: cnt['nonint'] += 1; continue
                vals = list(R.values())
                if len(set(vals)) != len(vals): cnt['dup'] += 1; continue
                if any(v % 4 for v in vals): cnt['mod4'] += 1; continue
                cc = coprime_class(vals)
                if cc is None: cnt['nocoprime'] += 1; ex = (a,b,c,lam,R); continue
                cnt['ok'] += 1
print(cnt, ex)
