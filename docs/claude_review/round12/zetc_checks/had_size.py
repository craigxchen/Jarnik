import numpy as np, itertools, sys
def sylvester(k):
    H = np.array([[1]])
    for _ in range(k): H = np.block([[H, H], [H, -H]])
    return H
def count(H, L):
    n = H.shape[0]; d = n-1
    # columns p = rows 1..n-1 of H (nonconstant); beta_p = H[p]; <k,beta> with k in Z^n_0, k_n = -sum k'
    B = H[1:, :]            # (n-1) x n
    Bp = B[:, :d] - B[:, d:d+1]   # b_{p,x} = beta_{p,x} - beta_{p,n}
    w = 2*L
    tot = 0; hist = {}
    # enumerate k' in {-1,0,1}^d in chunks over first 8 coords
    first = np.array(list(itertools.product((-1,0,1), repeat=8)))    # 6561 x 8
    rest = np.array(list(itertools.product((-1,0,1), repeat=d-8)))
    vf = first @ Bp[:, :8].T     # 6561 x (n-1)
    vr = rest @ Bp[:, 8:].T
    for i in range(len(rest)):
        v = vf + vr[i]
        ok = np.all(np.abs(v) <= w, axis=1)
        tot += int(ok.sum())
    return tot
for k, L in [(4, 1), (4, 2), (4, 3)]:
    H = sylvester(k)
    print('Sylvester n=%d L=%d |A|=%d' % (2**k, L, count(H, L)), flush=True)
