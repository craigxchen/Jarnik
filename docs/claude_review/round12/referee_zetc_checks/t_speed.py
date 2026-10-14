import time, numpy as np, flint
p = 1048571
n = 3000
X = np.random.randint(0, p, size=(n, n+100))
t = time.time(); M = flint.nmod_mat(X.tolist(), p); t1 = time.time()
r = M.rank(); t2 = time.time()
print("build", t1-t, "rank", r, t2-t1)
