import numpy as np, itertools, sys
from inc_cert import certify_inc
MOD = 0b10011
def mul(a, b):
    r = 0
    while b:
        if b & 1: r ^= a
        b >>= 1; a <<= 1
        if a & 16: a ^= MOD
    return r
pts = [1]
for i in range(14): pts.append(mul(pts[-1], 2))
def tr(y):
    y2 = mul(y, y); y4 = mul(y2, y2); y8 = mul(y4, y4); return y ^ y2 ^ y4 ^ y8
P = np.array([[tr(mul(a, x)) for a in range(1, 16)] for x in pts], dtype=np.int16)
tail = np.array(list(itertools.product((-1, 0, 1), repeat=10)), dtype=np.int8)
tailS = tail.astype(np.int16) @ P[5:]
L = 2; chunks = []
for head in itertools.product((-1, 0, 1), repeat=5):
    hs = np.array(head, dtype=np.int16) @ P[:5]
    ok = np.all(np.abs(tailS + hs) <= L, axis=1)
    if ok.any():
        sub = tail[ok]; chunks.append(np.concatenate([np.tile(np.array(head, dtype=np.int8), (len(sub), 1)), sub], 1))
A = np.concatenate(chunks); print("Syl16 L2 |A|", len(A), flush=True)
# Frobenius x -> x^2: point alpha^i -> alpha^(2i mod 15); check A invariant and phi sigma phi^-1 = sigma^2
frob = [(2 * i) % 15 for i in range(15)]           # point index i -> frob[i]
inv = np.argsort(frob)
Aphi = A[:, inv]                                     # (phi k)_{frob[i]} = k_i
S1 = set(map(bytes, A)); assert all(bytes(r) in S1 for r in Aphi), "A not Frobenius-invariant"
print("Frobenius invariance: OK", flush=True)
perm = [(i - 1) % 15 for i in range(15)]
reps = [(0,0),(0,1),(1,0),(1,1),(3,0),(3,1),(5,0),(5,1),(7,0),(7,1)]
rest = [(j, e) for e in (0, 1) for j in range(15) if (j, e) not in reps]
which = [tuple(map(int, a.split(","))) for a in sys.argv[1:]] if len(sys.argv) > 1 else reps + rest
certify_inc(A, perm, 15, 15, "Syl16L2_D15", which, batch=400)
