import pickle,sys
from cuttypes import CUTS
from gb import run,cuts_to_columns
T=pickle.load(open('types_7_4.pkl','rb'))
for j in [0,1]:
    cols=cuts_to_columns(T[j],CUTS)
    for drop in range(1,5):
        keep=[i for i in range(5) if i!=drop]
        c4=[(tuple(a[i] for i in keep)+(0,),e) for a,e in cols]
        if any(len(set(a[:4]))==1 for a,e in c4): print(j,drop,"constant column; skip"); continue
        out,tm=run(c4,4,32003,tag=f"C4_{j}_{drop}",timeout=900,slices=3,npts=4)
        print(j,drop,[ "".join(map(str,a[:4])) for a,e in c4],f"{tm:.1f}s",out.strip().replace("\n"," | ")[:200],flush=True)
