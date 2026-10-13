"""Compare referee scan with author's s9.pkl (pure data; opcodes checked) on N <= bound."""
import pickle, sys
from scan_ref_canon import canon
bound = int(float(sys.argv[2]))
mine = pickle.load(open(sys.argv[1], 'rb'))
auth = pickle.load(open('../L1L2_checks/s9.pkl', 'rb'))
A = {canon(cl): (N, C) for (N, fac, cl, C) in auth if N <= bound}
B = {canon(cl): (N, C) for (N, fac, cl, C) in mine if N <= bound}
print("author N<=bound:", len(A), " referee:", len(B))
onlyA = [A[k] for k in A if k not in B]; onlyB = [B[k] for k in B if k not in A]
print("only author:", len(onlyA), sorted(onlyA)[:10]); print("only referee:", len(onlyB), sorted(onlyB)[:10])
dC = max((abs(A[k][1] - B[k][1]) for k in A if k in B), default=0)
print("max |C diff| on common:", dC)
