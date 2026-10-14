import numpy as np, itertools
from paley_base import build, COLS
A=build(2)
V=np.load("paley_V.npy"); w=np.load("paley_w.npy")
# congruences mod 3
cong=[]
for c in itertools.product([0,1,2],repeat=11):
    c=np.array(c)
    if not c.any(): continue
    s=c@COLS
    if np.all(s%3==0): cong.append(c)
print("# nonzero congruence vectors mod 3:",len(cong), "(3^k-1 => k=%.3f)"%np.log(len(cong)+1)*1 if False else len(cong))
supp=[int((c>0).sum()) for c in cong]
print("support hist",np.bincount(supp).tolist())
ex=[c for c in cong if (c>0).sum()==5][:3]
for c in ex:
    s=c@COLS; print("c",c.tolist(),"sum",s.tolist())
    # lift with c in {-1,0,1}
    cl=np.where(c==2,-1,c); z=(cl@COLS)
    print("  lift -1/0/1 sum",z.tolist(),"div3 ok",np.all(z%3==0))
# width-7 directions, support 5
for s in [5,7,10]:
    sel=V[(w==7)&(np.abs(V).sum(1)==s)]
    print("width7 support",s,"example",sel[0].tolist(), "entries",sorted(set(sel.flatten().tolist())))
