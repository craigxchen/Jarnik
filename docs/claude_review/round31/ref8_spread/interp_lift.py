"""Reimplementation (referee's own code) of the deterministic interpolation completion behind Astra's
reported maximum 7,558, q = 383, S = first 189 positive rows of column 0.  tau' := the polynomial of degree
< |S| interpolating tau on S, evaluated on all of F_q (centred); then K^T tau' = 0 mod q by degree count."""
import numpy as np
q = 383
sq = {x*x % q for x in range(1, q)}
chi = np.array([0] + [1 if x in sq else -1 for x in range(1, q)])
K = chi[(np.arange(q)[None, :] - np.arange(q)[:, None]) % q]
A = [i for i in range(q) if K[i, 0] == 1]; S = A[:-2]
cen = lambda x: (x + q//2) % q - q//2
# Lagrange basis values l_i(j) mod q
L = np.zeros((len(S), q), dtype=object)
for a, i in enumerate(S):
    den = 1
    for k in S:
        if k != i: den = den * (i - k) % q
    dinv = pow(den, -1, q)
    for j in range(q):
        num = 1
        for k in S:
            if k != i: num = num * (j - k) % q
        L[a, j] = num * dinv % q
L = L.astype(np.int64)
mx = 0
for s in range(1, q//4 + 1):
    d = np.array([cen(pow(j, s, q)) for j in range(q)], dtype=np.int64)
    Kd = K @ d; assert (Kd % q == 0).all(); tau = Kd // q
    tp = np.array([cen(int(x)) for x in (tau[S] @ L) % q], dtype=np.int64)
    tp[S] = tau[S]
    assert ((K.T @ tp) % q == 0).all()
    lift = d - K.T @ tp
    assert ((lift - d) % q == 0).all() and len(set((K[S][:, 1:] @ lift[1:]).tolist())) == 1
    mx = max(mx, int(np.abs(lift).max()))
print("deterministic interpolation lifts s=1..95 OK; max |coefficient| =", mx)
