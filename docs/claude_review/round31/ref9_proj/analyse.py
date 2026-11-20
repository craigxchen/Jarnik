import json, sys
from math import log
from collections import defaultdict
for f in sys.argv[1:]:
    st = json.load(open(f))
    rows = st['rows']
    print("==", f, "windows", len(rows))
    byM = defaultdict(list)
    for r in rows: byM[r['M']].append(r)
    for M in sorted(byM):
        L = [r for r in byM[M] if r['Np'] > 1]
        if not L: continue
        n = len(L)
        t1 = sum(r['minpairT'] == 1 for r in L)
        tP = sum(r['minpairT'] is not None and r['minpairT'] <= r['Pp']**2 for r in L)
        fr = sorted(log(r['minT'])/log(r['Np']) if r['minT'] > 1 else 0.0 for r in L)
        frp = sorted(log(r['minpairT'])/log(r['Np']) if r['minpairT'] > 1 else 0.0 for r in L)
        med = fr[n//2]; medp = frp[n//2]
        # certificate beats the trivial projection bound? (1+l1)^2 T < M^2 means X_b small; record min ratio
        mc = min(r['minCert'] for r in L); mp = min(r['bestproj'] for r in L)
        print(f"M={M:2d} n={n:6d}  pairT=1: {t1/n:6.3f}  pairT<=P+^2: {tP/n:6.3f}  median log minT/log N' (all c): {med:5.3f} (pairs: {medp:5.3f})  "
              f"min (1+|c|)sqrtT/M^1.5={mc:6.3f}  min (2+X_b)/M^1.5={mp:6.3f}")
