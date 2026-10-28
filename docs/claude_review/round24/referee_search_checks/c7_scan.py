import mpmath as mp, itertools, sys
mp.mp.dps=40
from c4_construct import curve, places, ray
import numpy as np
N=int(sys.argv[1]); worst=[]
bad=0
for i in range(1,N):
    m=-3+2*mp.mpf(i)/N
    try:
        C=curve(m); pl=places(C)
        rays=sorted(ray(v) for w,v in pl)
        ok=(len(pl)==30 and len(set(rays))==30)
        ws=[w for w,v in pl]
        sep=min(abs(a-b) for a,b in itertools.combinations(ws,2))/max(abs(w) for w in ws)
        rd=min(abs(mp.im(w))/abs(w) for w in ws)
        rec=(float(m),ok,float(sep),float(rd),float(mp.re(C['e'])),float(C['D']),float(mp.re(C['e'])-C['D']),float(C['sig']))
    except Exception as ex:
        rec=(float(m),False,str(ex)); bad+=1
    worst.append(rec)
    if not rec[1]: print('PROBLEM',rec)
oks=[r for r in worst if r[1]]
print('N',N,'all ok:',len(oks)==len(worst),' min rel sep',min(r[2] for r in oks),' min rel dist to R',min(r[3] for r in oks))
print('sign ranges e:',min(r[4] for r in oks),max(r[4] for r in oks),' D:',min(r[5] for r in oks),max(r[5] for r in oks),' e-D:',min(r[6] for r in oks),max(r[6] for r in oks),' sig:',min(r[7] for r in oks),max(r[7] for r in oks))
# rel sep as function of m near ends
for r in oks[:3]+oks[-3:]: print(r)
