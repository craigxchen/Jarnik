import pickle,sys
from cuttypes import CUTS
from gb import run,cuts_to_columns
D,h,char=map(int,sys.argv[1:4])
js=sys.argv[4] if len(sys.argv)>4 else "all"
T=pickle.load(open(f'types_{D}_{h}.pkl','rb'))
rng=range(len(T)) if js=="all" else [int(x) for x in js.split(",")]
for j in rng:
    cols=cuts_to_columns(T[j],CUTS)
    out,tm=run(cols,h,char,tag=f"S{D}_{h}_{j}_{char}",timeout=int(sys.argv[5]) if len(sys.argv)>5 else 3600)
    print(D,h,j,["".join(map(str,a)) for a,e in cols],f"{tm:.1f}s",out.strip().replace("\n"," | ")[:300],flush=True)
