import pickle, sympy
from cuttypes import CUTS
from gb import cuts_to_columns
T=pickle.load(open('sqfree_types.pkl','rb'))
for key,L in T.items():
    rk=[]
    for ty in L:
        cols=cuts_to_columns(ty,CUTS)
        M=sympy.Matrix([[a[i]-a[0] for a,e in cols] for i in range(1,5)])
        rk.append(M.rank())
    print(key,len(L),sorted(set(rk)))
