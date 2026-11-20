"""Independent verification of an integral-cofactor configuration as posed in
reply11.txt (A2), from first principles: rows, H_ij, Q_ij, t_ij.
Usage: python3 verify.py P C i j k [...]   (indices into the eta=1 mark list)."""
import math, pickle, sys
from marks import gmul

P = int(sys.argv[1]); C = int(sys.argv[2]); idx = [int(a) for a in sys.argv[3:]]
ps, pis, q, marks = pickle.load(open(f"marks_{P}.pkl", "rb"))
s = len(ps)

def conj(z): return (z[0], -z[1])
def f(r, k): return pis[k] if (r >> k) & 1 else conj(pis[k])

# exact checks of the canonical identity A_j - A_i = 2 i H T, T = -Im X
A = {r: marks[r] for r in idx}
G = 0
data = {}
for a in range(len(idx)):
    for b in range(a + 1, len(idx)):
        i, j = idx[a], idx[b]
        H = (1, 0); X = (1, 0); Q = 1
        for k in range(s):
            if ((i >> k) & 1) == ((j >> k) & 1):
                H = gmul(H, f(i, k))
            else:
                X = gmul(X, f(i, k)); Q *= ps[k]
        T = -X[1]
        lhs = (A[j][0] - A[i][0], A[j][1] - A[i][1])
        rhs = gmul((0, 2), gmul(H, (T, 0)))
        assert lhs == rhs, "canonical identity"
        assert 0 < T * T < Q
        assert H[0] ** 2 + H[1] ** 2 == q // Q
        data[(i, j)] = (H, Q, T)
        G = math.gcd(G, T)
print("canonical identities exact; G = gcd T_ij =", G)
lam_den = G  # lambda = 1/G
w_num = (-A[idx[0]][0], -A[idx[0]][1])  # w = -A_1/G, so z_1 = 0
z = {}
for r in idx:
    nx, ny = A[r][0] + w_num[0], A[r][1] + w_num[1]
    assert nx % G == 0 and ny % G == 0, "z not integral"
    z[r] = (nx // G, ny // G)
print("z =", z)
ok = True
for (i, j), (H, Q, T) in data.items():
    t = T // G
    assert t * G == T and t != 0
    diff = (z[j][0] - z[i][0], z[j][1] - z[i][1])
    assert diff == gmul((0, 2), gmul(H, (t, 0))), "fake difference form"
    # |t| <= (C/2) sqrt(Q) q^{-1/4}  <=>  16 t^4 q <= C^4 Q^2
    good = 16 * t ** 4 * q <= C ** 4 * Q * Q
    ok &= good
    print(f"pair {i},{j}: t={t}, log2 Q={math.log2(Q):.2f}, |t|/((1/2)sqrt(Q)q^-1/4)={abs(t)/(0.5*math.sqrt(Q)*q**-0.25):.4f}, ok={good}")
    print("   norms |z|^2:", z[i][0]**2 + z[i][1]**2, z[j][0]**2 + z[j][1]**2)
# triangle identities are automatic: z are actual Gaussian integers
print("ALL CONDITIONS (A2) at C =", C, ":", ok)
