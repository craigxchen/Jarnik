import json, sys
from spec import *
for f in sys.argv[1:]:
    d=json.load(open(f)); N=d['bezout_N']
    norms={}; bad=[]; cas=[]; asc=[]; allNdiv=True; maxg=0; oddR=True; unit1=0; pmall=0
    for T in [x for x in range(-300,3001) if abs(x)>1]:
        zs=[geval(c,T) for c in d['coeffs']]
        if any(z==(0,0) for z in zs): bad.append(('zero',T)); continue
        r=analyse(d,T)
        norms[r['gcd_norm']]=norms.get(r['gcd_norm'],0)+1
        allNdiv &= r['g_divides_N']
        assert r['same_norm'] and r['primitive']
        if not r['distinct']: bad.append(T)
        if r['associates']: asc.append(T)
        if r['conj_associates']: cas.append(T)
        if r['units'] and len(set(r['units']))==1: unit1+=1
        if r['pm_pairs_all']: pmall+=1
    print(d['n'], d['seeds'], 'N=',N, 'all g|N:',allNdiv, 'gcd norms:',dict(sorted(norms.items())), 'non-distinct T:',bad, 'assoc T:',asc,'conj-assoc T:',cas, 'one-unit-class(primary) count:',unit1,'all pairs +-P/Pbar count:',pmall)
