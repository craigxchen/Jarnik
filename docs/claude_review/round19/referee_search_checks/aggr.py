import json, sys
from collections import Counter
f=sys.argv[1]
Ds=sorted(set(json.loads(l)['D'] for l in open(f)))
for D in Ds:
    hg,hp,hs=Counter(),Counter(),Counter(); cnt=Counter(); mism=[]; sp=[]; cpu=0; top=[]; chunks=0
    for l in open(f):
        r=json.loads(l)
        if r['D']!=D: continue
        chunks+=1
        for k,v in r['hist_generic'].items(): hg[int(k)]+=v
        for k,v in r['hist_pm'].items(): hp[int(k)]+=v
        for k,v in r['hist_special'].items(): hs[int(k)]+=v
        cnt.update(r['counts']); mism+=r['mismatches']; sp+=r['specials']; cpu+=r['cpu']; top+=r['top']
    print('D',D,'chunks',chunks, dict(cnt), 'cpu',round(cpu,1))
    print(' pm', sorted(hp.items(),reverse=True)[:8], 'total', sum(hp.values()), 'distinct', len(hp))
    print(' gen', sorted(hg.items(),reverse=True)[:8], 'total', sum(hg.values()))
    print(' spec hist', dict(hs), 'n', len(sp), 'mism', len(mism), mism[:3])
    print(' special lambdas', Counter(str(s[3]) for s in sp))
    print(' top>-40 count', Counter(t[0] for t in top))
