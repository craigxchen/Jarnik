import numpy as np, pickle, itertools
from newton import search
from cuttypes import CUTS
from gb import cuts_to_columns
T=pickle.load(open('types_9_5.pkl','rb'))
cols=cuts_to_columns(T[0],CUTS)
found,deg,fail=search(cols,5,3000,seed=7)
print(len(found),deg,fail)
pickle.dump((cols,[z for d,z in found]),open('d9_found.pkl','wb'))
n=len(cols)
for d,z in found[:40]:
    r,s=z[:n],z[n:]
    # reality test under gauge r0-s0=1: real configs have s=-conj(r) after translation; test via the invariant: multiset {s_k} vs {-conj r_k} up to translation
    shift=np.mean(r+s)  # zero by gauge
    real_err=np.max(np.abs(s+np.conj(r)))
    print(f"sep {d:.4f} real_err {real_err:.3e}  r={np.round(r,4)}")
