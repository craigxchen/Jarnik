from comb import *
import io,contextlib
def vec(pos,neg):
    v=[0]*8
    for j in pos: v[j]=1
    for j in neg: v[j]=-1
    return tuple(v)
SEL1=[vec((0,4),(1,3)),vec((0,6),(1,7)),vec((1,5),(3,4)),vec((2,6),(3,4)),vec((3,6),(4,7))]
SEL2=[vec((0,3),(1,4)),vec((0,3),(2,7)),vec((1,4),(2,7)),vec((0,4),(1,7)),vec((1,3),(2,4)),vec((3,6),(4,5)),vec((4,6),(5,7))]
for name,roots,sel in (('C1',C1,SEL1),('C2',C2,SEL2)):
    with contextlib.redirect_stdout(io.StringIO()):
        dirs=analyse(name,roots)
    idx=[]
    for s in sel:
        m=[k for k,d in enumerate(dirs) if d==s or d==tuple(-a for a in s)]
        assert len(m)==1; idx.append(m[0])
    print(name,'selected indices in my order',idx, 'b values',[Bpoly(roots,dirs[k])[0][1] for k in idx])
    analyse(name+' selected',roots,idx)
