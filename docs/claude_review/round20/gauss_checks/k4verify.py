"""Exact verification (Python ints) of the k=4 degree-1 family H_T(s) = s - r_T,
 r_123 = i, r_12 = -6-4i, r_13 = -6-3i, r_23 = -6+5i, r_1 = -7+6i, r_2 = -3-2i, r_3 = -4-3i.
Polynomial identities: Im P_1 = -420, Im P_2 = -180, Im P_3 = -240 (checked symbolically in k4family/k4verify_sym).
For integer s: lattice points z_x (x=0..3), common gcd g (bounded: only primes dividing resultants), reduced radius R',
reduced max chord c', ratio c'/R'^(3/7) and exponent log c'/log R'; anchored residues Y_i, t_ij = det(P_i,P_j)/G_ij with
G_ij = Norm gcd(P_i,P_j) computed directly (Gaussian gcd), and the block-product prediction prod_{T∋i,j} n_T."""
import math
from itertools import combinations
from gcdred import ggcd, gnorm, gmul, gexdiv
def conj(a): return (a[0], -a[1])
roots = {(1,2,3): (0,1), (1,2): (-6,-4), (1,3): (-6,-3), (2,3): (-6,5), (1,): (-7,6), (2,): (-3,-2), (3,): (-4,-3)}
def data(sv):
    h = {T: (sv - r[0], -r[1]) for T, r in roots.items()}
    z = []
    for x in (0,1,2,3):
        acc = (1,0)
        for T in h: acc = gmul(acc, h[T] if x in T else conj(h[T]))
        z.append(acc)
    N = gnorm(z[0]); assert all(gnorm(w) == N for w in z)
    g = z[0]
    for w in z[1:]: g = ggcd(g, w)
    Zr = [gexdiv(w, g) for w in z]; Np = N // gnorm(g)
    assert all(gnorm(w) == Np for w in Zr)
    P = {}
    for i in (1,2,3):
        acc = (1,0)
        for T in h:
            if i in T: acc = gmul(acc, h[T])
        P[i] = acc
    Y = [P[i][1] for i in (1,2,3)]
    G = {(i,j): gnorm(ggcd(P[i], P[j])) for i,j in combinations((1,2,3),2)}
    pred = {(i,j): math.prod(gnorm(h[T]) for T in h if i in T and j in T) for i,j in combinations((1,2,3),2)}
    t = {(i,j): (P[i][0]*P[j][1]-P[j][0]*P[i][1]) // G[(i,j)] for (i,j) in G}
    assert all((P[i][0]*P[j][1]-P[j][0]*P[i][1]) % G[(i,j)] == 0 for (i,j) in G)
    ch2 = max((Zr[a][0]-Zr[b][0])**2+(Zr[a][1]-Zr[b][1])**2 for a, b in combinations(range(4), 2))
    lch = 0.5*math.log(ch2); lR = 0.5*math.log(Np)
    blocks = [gnorm(h[T]) for T in h]
    eta = max(abs(math.log(n) - 2*math.log(sv)) for n in blocks) / (2*math.log(sv))
    return dict(s=sv, Ngcd=gnorm(g), Y=Y, t=list(t.values()), Gratio=[G[k]//pred[k] if G[k]%pred[k]==0 else ('!',G[k],pred[k]) for k in G],
                ratio=math.exp(lch - (3/7)*lR), expo=lch/lR, eta=eta)
if __name__ == "__main__":
    for sv in [101, 1001, 10**4+1, 10**5+3, 10**6+1, 10**9+7, 10**15+37, 10**30+57]:
        d = data(sv)
        print("s=%s Norm(g)=%d Y=%s t=%s G/prod n_T=%s c'/R'^(3/7)=%.5f exponent=%.6f eta=%.2e" % (
            d['s'], d['Ngcd'], d['Y'], d['t'], d['Gratio'], d['ratio'], d['expo'], d['eta']))
