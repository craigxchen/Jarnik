import json,collections,sys
for fn in sys.argv[1:]:
    hp=collections.Counter(); hg=collections.Counter(); hs=collections.Counter(); st=collections.Counter(); sp=[]; n=0; tm=0; top=[]; keys=set()
    for l in open(fn):
        r=json.loads(l); k=(r['D'],r['orbit'],r['c0'])
        if k in keys: continue
        keys.add(k); n+=r['n']; tm+=r['time']
        hp.update({int(k):v for k,v in r['hist_pm'].items()}); hg.update({int(k):v for k,v in r['hist_generic'].items()}); hs.update({int(k):v for k,v in r['hist_special'].items()}); st.update(r['stats']); sp+=r['special']; top+=r['top']
    print(fn,'chunks',len(keys),'pairs',n,'cpu time',round(tm,1))
    print('  lam=+-1 hist top',sorted(hp.items(),reverse=True)[:8],'total',sum(hp.values()))
    print('  generic hist top',sorted(hg.items(),reverse=True)[:5])
    print('  special (non-sqfree, lam!=+-1) hist',sorted(hs.items(),reverse=True),'count',len(sp))
    print('  stats',dict(st)); print('  top>-40 count',len(top),top[:4])
