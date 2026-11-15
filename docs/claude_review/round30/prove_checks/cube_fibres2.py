# Correct hashing (round to 1e-30). Max fibre (= max cluster over all base points) at generic members and at the
# special real points found by cube_fibres.py (parsed from its output).
import itertools, mpmath as mp, random, re, sys
from collections import Counter
mp.mp.dps=60
eps=[1,-1,-1,1]
def fibres(yv):
    roots=[]
    for k in range(4):
        r=mp.mpc(eps[k]/yv[k],yv[k]); roots+=[r,-r]
    C=Counter(); members={}
    for s in itertools.product([1,-1],repeat=8):
        v=tuple(int(mp.nint(sum(s[j]*mp.im(roots[j]**m) for j in range(8))*mp.mpf(10)**30)) for m in (1,2,3))
        C[v]+=1
    return sorted(C.values())[-4:]
random.seed(2)
for trial in range(4):
    while True:
        a=mp.mpf(random.uniform(-3,3)); b=mp.mpf(random.uniform(-3,3))
        disc=(a*b*(a+b))**2-4*a*b
        if disc>0: break
    c=(-a*b*(a+b)+mp.sqrt(disc))/(2*a*b); yv=[a,b,c,-(a+b+c)]
    print('generic',[mp.nstr(x,6) for x in yv],fibres(yv))
txt=open(sys.argv[1]).read()
# refine special points from printed 8 digits by Newton on the system (sum d y=0, sum d g=0, prod=1)
pts=re.findall(r"\(([-\d, ]+)\) \['([-\d.]+)', '([-\d.]+)', '([-\d.]+)', '([-\d.]+)'\]",txt)
done=set()
for p in pts:
    d=[int(x) for x in p[0].split(',')]; y0=[mp.mpf(x) for x in p[1:]]
    f=lambda a,b,c: [sum(d[k]*v for k,v in enumerate([a,b,c,-(a+b+c)])),
                     sum(d[k]*(3/v-v**3) for k,v in enumerate([a,b,c,-(a+b+c)])),
                     a*b*c*(-(a+b+c))-1]
    sol=mp.findroot(f,y0[:3])
    yv=[sol[0],sol[1],sol[2],-(sol[0]+sol[1]+sol[2])]
    key=tuple(mp.nstr(v,10) for v in yv)
    if key in done: continue
    done.add(key)
    print(d,[mp.nstr(v,12) for v in yv],'fibres',fibres(yv))
