import pickle,sys
from gb import run
D,h,char=map(int,sys.argv[1:4])
T=pickle.load(open(f'nsf_{D}_{h}.pkl','rb'))
for j,cols in enumerate(T):
    out,tm=run(list(cols),h,char,tag=f"N{D}_{h}_{j}_{char}",timeout=1800)
    print(D,h,j,[("".join(map(str,a)),e) for a,e in cols],f"{tm:.1f}s",out.strip().replace("\n"," | ")[:300],flush=True)
