import json, sys
from collections import Counter
D=int(sys.argv[1]); f=sys.argv[2]
AU='/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round19/search_checks/'
d=json.load(open(AU+'out/D%d_analysis.json'%D))
def todict(h):
    if isinstance(h,dict): return {int(k):v for k,v in h.items()}
    return {int(k):v for k,v in h}
hp=Counter(); hg=Counter(); me=set()
for l in open(f):
    r=json.loads(l)
    for k,v in r['hist_pm'].items(): hp[int(k)]+=v
    for k,v in r['hist_generic'].items(): hg[int(k)]+=v
    for s in r['specials']: me.add((frozenset(map(tuple,s[1])),tuple(sorted(map(tuple,s[2]))),tuple(s[3])))
au_sp=set((frozenset(map(tuple,r['A'])),tuple(sorted(map(tuple,r['C']))),tuple(r['lam'])) for r in d['records_gain_gt_-40_and_all_special'] if tuple(r['lam']) not in ((1,1),(-1,1)))
print('D',D,'pm hist equal',todict(d['hist_lam_pm1'])==dict(hp),'distinct',len(hp),'| generic equal',todict(d['hist_generic'])==dict(hg),'| specials author',len(au_sp),'mine',len(me),'same set',au_sp==me)
print(' example author hist head', str(d['hist_lam_pm1'])[:120])
