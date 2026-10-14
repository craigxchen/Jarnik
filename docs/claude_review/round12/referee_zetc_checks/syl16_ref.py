import numpy as np, itertools, sys, time
from cyc_rank import certify
MOD = 0b10011
def mul(a, b):
    r = 0
    while b:
        if b & 1: r ^= a
        b >>= 1; a <<= 1
        if a & 16: a ^= MOD
    return r
alpha = 2
pts = [1]
for i in range(14): pts.append(mul(pts[-1], alpha))
assert len(set(pts)) == 15
def tr(y):
    y2 = mul(y, y); y4 = mul(y2, y2); y8 = mul(y4, y4)
    t = y ^ y2 ^ y4 ^ y8; assert t in (0, 1); return t
P = np.array([[tr(mul(a, x)) for a in range(1, 16)] for x in pts], dtype=np.int16)  # 15 points x 15 columns
assert all(P[:, c].sum() == 8 for c in range(15))
d = 15
tail = np.array(list(itertools.product((-1, 0, 1), repeat=10)), dtype=np.int8)
tailS = tail.astype(np.int16) @ P[5:]
L_list = [int(a) for a in sys.argv[1:]] or [1, 2, 3]
for L in L_list:
    chunks = []
    for head in itertools.product((-1, 0, 1), repeat=5):
        hs = np.array(head, dtype=np.int16) @ P[:5]
        ok = np.all(np.abs(tailS + hs) <= L, axis=1)
        if ok.any():
            sub = tail[ok]
            chunks.append(np.concatenate([np.tile(np.array(head, dtype=np.int8), (len(sub), 1)), sub], 1))
    A = np.concatenate(chunks)
    print("Syl16 L", L, "|A|", len(A), flush=True)
    if L == 2 and "--cert" in sys.argv[0:0] or (L == 2 and len(sys.argv) > 1 and sys.argv[1] == "2"):
        perm = [(i - 1) % 15 for i in range(15)]
        reps = [(0,0),(0,1),(1,0),(1,1),(3,0),(3,1),(5,0),(5,1),(7,0),(7,1)]
        rest = [(j, e) for e in (0, 1) for j in range(15) if (j, e) not in reps]
        certify(A, perm, 15, 15, "Syl16L2_D15", chis=reps + rest)
