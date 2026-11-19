import pickle,sys
from cuttypes import CUTS
from gb import cuts_to_columns
from ms import run_sing
D,h,char,j,sl=map(int,sys.argv[1:6])
T=pickle.load(open(f'types_{D}_{h}.pkl','rb')) if (D,h)!=(0,0) else None
cols=cuts_to_columns(T[j],CUTS)
out,tm=run_sing(cols,h,char,f"o2_{D}_{h}_{j}_{char}_{sl}",slices=sl,timeout=1500)
print(D,h,j,["".join(map(str,a)) for a,e in cols],f"{tm:.1f}s",out,flush=True)
