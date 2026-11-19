import pickle,sys
from cuttypes import enum,CUTS
from gb import run,cuts_to_columns
D,h,char=int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3])
T=enum(D,h)
print(len(T),"types")
pickle.dump(T,open(f'types_{D}_{h}.pkl','wb'))
lim=int(sys.argv[4]) if len(sys.argv)>4 else 10**9
for j,ty in enumerate(T[:lim]):
    cols=cuts_to_columns(ty,CUTS)
    desc=[ "".join(map(str,a)) for a,e in cols]
    out,tm=run(cols,h,char,tag=f"c{D}_{h}_{j}_{char}",timeout=600,slices=int(sys.argv[5]))
    print(D,h,j,desc,f"{tm:.1f}s",out.strip().replace("\n"," | ")[:300],flush=True)
