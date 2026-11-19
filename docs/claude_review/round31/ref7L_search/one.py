import pickle,sys
from cuttypes import CUTS
from gb import run,cuts_to_columns
D,h,char,j,sl=map(int,sys.argv[1:6])
T=pickle.load(open(f'types_{D}_{h}.pkl','rb'))
cols=cuts_to_columns(T[j],CUTS)
out,tm=run(cols,h,char,tag=f"o{D}_{h}_{j}_{char}",timeout=1200,slices=sl)
print(D,h,j,["".join(map(str,a)) for a,e in cols],f"{tm:.1f}s",out.strip().replace("\n"," | ")[:400],flush=True)
