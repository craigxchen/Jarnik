import time, math, collections
from m8lib import *
t=time.time()
A, cols = m8_box(2)
print("|A| =", len(A))
H, JA = hilbert_JA(A, 8, exact=True)
print("H_A (exact Q) =", H, "sum", sum(H)); print("dim J(A)_i (exact Q) =", JA, time.time()-t)
W = width_set(A, R=3, maxdeg=8)
print("#primitive directions |v_i|<=3 with tau+1<=8:", len(W))
print("count by (tau,h):", sorted(collections.Counter((t,h) for v,t,h in W).items()))
print("any tau<h:", [x for x in W if x[1]<x[2]][:10])
print("max |v_i| among them:", max(max(abs(z) for z in v) for v,_,_ in W))
import pickle; pickle.dump((A,cols,H,JA,W), open('m8_step1.pkl','wb'))
