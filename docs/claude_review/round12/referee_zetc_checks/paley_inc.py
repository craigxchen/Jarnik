import numpy as np, itertools
from inc_cert import certify_inc
q = 11; D0 = [0,1,3,4,5,9]
cols = [sorted((i + x) % q for x in D0) for i in range(q)]
P = np.array([[1 if x in c else 0 for c in cols] for x in range(q)])
allk = np.array(list(itertools.product((-1,0,1), repeat=11)), dtype=np.int8)
A = allk[np.all(np.abs(allk.astype(int) @ P) <= 2, axis=1)]
print("|A|", len(A))
perm = [(i - 1) % q for i in range(q)]
certify_inc(A, perm, 8, q, "PaleyL2_D8_inc", [(j, e) for e in (0, 1) for j in range(q)])
certify_inc(A, perm, 7, q, "PaleyL2_D7_inc", [(0, 0)])
