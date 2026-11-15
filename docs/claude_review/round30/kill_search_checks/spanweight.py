# kill_search: for each code class, min l1-weight of integer vectors in the real span (coeffs in [-1,1] over the code's rows,
# plus a brute-force search over {-1,0,1}^D). Lemma 1.1 needs weight >= L+1.
import sys, pickle, itertools, numpy as np
fn, D = sys.argv[1], int(sys.argv[2]); L = D-1-D//2
codes = pickle.load(open(fn, 'rb'))
W = np.array(list(itertools.product([-1,0,1], repeat=D)), np.int8)
for key in codes:
    P = np.array(key).T.astype(float)
    B = P[1:] - P[0]
    rk = np.linalg.matrix_rank(B)
    U, S, Vt = np.linalg.svd(B); Q = Vt[:rk].T
    res = np.linalg.norm(W - (W @ Q) @ Q.T, axis=1)
    In = W[(res < 1e-9) & (np.abs(W).sum(1) > 0)]
    w = np.abs(In).sum(1)
    a = In[np.argmin(w)]
    print("rank", rk, "min weight of {0,+-1} vector in span:", w.min(), "(need >=", L+1, ")", "example", a.tolist(),
          "-> KILLED by Lemma 1.1" if w.min() < L+1 else "")
