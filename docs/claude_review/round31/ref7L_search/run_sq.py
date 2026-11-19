import pickle,sys
from cuttypes import CUTS
from gb import run,cuts_to_columns
T=pickle.load(open('sqfree_types.pkl','rb'))
D,h,char=int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3])
for j,ty in enumerate(T[(D,h)]):
    cols=cuts_to_columns(ty,CUTS)
    desc=[ "".join(map(str,a)) for a,e in cols]
    out,tm=run(cols,h,char,tag=f"{D}_{h}_{j}_{char}")
    print(D,h,j,desc,f"{tm:.1f}s",out.strip().replace("\n"," | ")[:600],flush=True)
