import numpy as np, itertools
def paley12():
    q = 11
    sq = {(x*x) % q for x in range(1, q)}
    chi = lambda a: 0 if a % q == 0 else (1 if a % q in sq else -1)
    # Paley I: H = I + S where S skew conference matrix of order 12
    n = 12
    S = np.zeros((n, n), dtype=int)
    for i in range(1, n):
        S[0, i] = 1; S[i, 0] = -1
    for i in range(q):
        for j in range(q):
            S[i+1, j+1] = chi(j - i)
    H = np.eye(n, dtype=int) + S
    assert (H @ H.T == n*np.eye(n)).all()
    # normalise so that column 0 (base point) ... make first row all ones: multiply columns by sign of row 0
    H = H * H[0:1, :]
    # make first column all ones? not needed: profile = rows 1..n-1 as sign vectors on points (columns)
    return H
def count(H, L):
    n = H.shape[0]; d = n-1
    B = H[1:, :]
    Bp = B[:, :d] - B[:, d:d+1]
    w = 2*L
    split = min(8, d)
    first = np.array(list(itertools.product((-1,0,1), repeat=split)))
    rest = np.array(list(itertools.product((-1,0,1), repeat=d-split)))
    vf = first @ Bp[:, :split].T
    vr = rest @ Bp[:, split:].T
    tot = 0
    for i in range(len(rest)):
        ok = np.all(np.abs(vf + vr[i]) <= w, axis=1); tot += int(ok.sum())
    return tot
H = paley12()
print('row sums (should be 0 for rows 1..):', H[1:].sum(axis=1))
for L in (1, 2, 3):
    print('Paley n=12 L=%d |A|=%d  (cube 3^11=%d)' % (L, count(H, L), 3**11), flush=True)
