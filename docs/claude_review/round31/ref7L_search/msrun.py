import pickle,sys,subprocess,time
from cuttypes import CUTS
from gb import cuts_to_columns
from ms import system
D,h,char,j,sl=map(int,sys.argv[1:6])
T=pickle.load(open(f'types_{D}_{h}.pkl','rb'))
cols=cuts_to_columns(T[j],CUTS)
v,e=system(cols,h,sl,1)
fn=f"m_{D}_{h}_{j}_{char}_{sl}"
open(fn+".in","w").write(",".join(v)+"\n"+str(char)+"\n"+",\n".join(e)+"\n")
t=time.time()
r=subprocess.run(["prlimit","--as=3500000000","msolve","-f",fn+".in","-o",fn+".out","-g","2","-t","1"],capture_output=True,text=True,timeout=3000)
print(r.stdout[-2000:],r.stderr[-2000:],time.time()-t)
