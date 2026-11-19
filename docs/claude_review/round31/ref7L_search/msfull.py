# msolve with per-pair Rabinowitsch variables killing coincidence components; prime field parametrisation
import pickle,sys,random,itertools,subprocess,time,sympy
from cuttypes import CUTS
from gb import cuts_to_columns
from ms import system
D,h,p,j=map(int,sys.argv[1:5])
T=pickle.load(open(f'types_{D}_{h}.pkl','rb'))
cols=cuts_to_columns(T[j],CUTS)
v,e=system(cols,h,0,1)
n=len(cols); rnd=random.Random(5)
T_eq=e[-1]; e=e[:-1]; v=v[:-1]   # drop T, re-add at end
extra=[]
for k,l in itertools.combinations(range(n),2):
    a,b=rnd.randint(2,9999),rnd.randint(2,9999)
    w=f"w{k}_{l}"; extra.append(w)
    ex=sympy.expand(sympy.sympify(f"{w}*((x{k}-x{l})+{a}*(y{k}+y{l}))*((x{k}-x{l})+{b}*(y{k}-y{l}))-1"))
    e.append(str(ex).replace("**","^").replace(" ",""))
vars_=v+extra+["T"]
e.append("T-"+"-".join(f"{rnd.randint(1,997)}*{x}" for x in v))
fn=f"full_{D}_{h}_{j}_{p}"
open(fn+".in","w").write(",".join(vars_)+"\n"+str(p)+"\n"+",\n".join(e)+"\n")
t=time.time()
r=subprocess.run(["prlimit","--as=3500000000","msolve","-f",fn+".in","-o",fn+".out","-t","1"],capture_output=True,text=True)
print("time",time.time()-t, r.stderr[-300:])
print(open(fn+".out").read()[:300])
