# many complex Newton starts for D=9 type; collect distinct nondegenerate solutions (rounded), save
import numpy as np, pickle, sys
from newton import search
from cuttypes import CUTS
from gb import cuts_to_columns
T=pickle.load(open('types_9_5.pkl','rb'))
cols=cuts_to_columns(T[0],CUTS)
allz=[]
for seed in range(int(sys.argv[1]),int(sys.argv[2])):
    found,deg,fail=search(cols,5,1000,seed=1000+seed)
    allz+=[z for d,z in found]
    keys=set(tuple(np.round(np.concatenate([z.real,z.imag]),6)) for z in allz)
    print(seed,len(found),"distinct so far",len(keys),flush=True)
pickle.dump((cols,allz),open(f'd9_all_{sys.argv[1]}.pkl','wb'))
