import numpy as np, itertools
from cyc_rank import certify
q = 11
QR = {(x*x) % q for x in range(1, q)}
chi = lambda a: 0 if a % q == 0 else (1 if a % q in QR else -1)
# Paley I (q=11, q=3 mod 4): H = I + S, S = [[0, 1^T], [-1, Q]], Q_ij = chi(j - i)
n = 12
S = np.zeros((n, n), int)
S[0, 1:] = 1; S[1:, 0] = -1
for i in range(q):
    for j in range(q):
        S[1+i, 1+j] = chi(j - i)
H = np.eye(n, dtype=int) + S
assert np.all(H @ H.T == n*np.eye(n))
# points = columns; base point = column 0 (infinity); normalise a row to all ones: row 0 already all ones? 
# make row 0 all ones by negating columns
sg = H[0].copy(); H = H * sg[None, :]
assert np.all(H[0] == 1) and np.all(H @ H.T == n*np.eye(n))
# patterns relative to base point column 0: S_p = {x>=1 : H[p,x] != H[p,0]}
pats = [np.array([1 if H[p, x] != H[p, 0] else 0 for x in range(1, n)]) for p in range(1, n)]
print("pattern sizes", [int(s.sum()) for s in pats])
d = n - 1
allk = np.array(list(itertools.product((-1, 0, 1), repeat=d)), dtype=np.int8)
P = np.array(pats).T   # d x 11
for L in (1, 2):
    A = allk[np.all(np.abs(allk.astype(int) @ P) <= L, axis=1)]
    print("Paley12 L", L, "|A|", len(A))
    # translation x -> x+1 on F_11: point index 1+i -> 1+((i+1)%11); (sigma k)_i = k_{perm[i]}
    perm = [((i - 1) % q) for i in range(q)]
    if L == 2:
        certify(A, perm, 8, q, "Paley12L2_D8")
        certify(A, perm, 7, q, "Paley12L2_D7", chis=[(0,0),(0,1),(1,0),(1,1)])
