import itertools, numpy as np, flint
def exact_H(A, maxdeg):
    d = len(A[0])
    mons = [e for t in range(maxdeg+1) for e in itertools.product((0,1,2), repeat=d) if sum(e) == t] if d <= 8 else None
    if mons is None:
        mons = []
        for t in range(maxdeg+1):
            for twos in range(t//2+1):
                ones = t - 2*twos
                for P2 in itertools.combinations(range(d), twos):
                    rest = [x for x in range(d) if x not in P2]
                    for P1 in itertools.combinations(rest, ones):
                        e = [0]*d
                        for x in P2: e[x] = 2
                        for x in P1: e[x] = 1
                        mons.append(tuple(e))
    H = []
    for t in range(maxdeg+1):
        cols = [e for e in mons if sum(e) <= t]
        M = flint.fmpz_mat([[int(np.prod([k[i]**e[i] for i in range(d) if e[i]])) if any(e) else 1 for e in cols] for k in A])
        H.append(M.rank())
    return H
def syl(kk):
    n = 2**kk
    pts = list(range(1, n))
    P = np.array([[bin(a & x).count("1") % 2 for a in range(1, n)] for x in pts])
    return P
def model(P, L):
    d = P.shape[0]
    allk = np.array(list(itertools.product((-1,0,1), repeat=d)), dtype=np.int8)
    return [tuple(map(int, k)) for k in allk[np.all(np.abs(allk.astype(int) @ P) <= L, axis=1)]]
# H_8 (Sylvester), L=1
P8 = syl(3)
for L in (1,):
    A = model(P8, L); print("H8 L", L, "|A|", len(A), "H", exact_H(A, 4))
for L in (2, 3):
    print("H8 L", L, "|A|", len(model(P8, L)))
# Paley 12 L=1
q = 11; QR = {(x*x) % q for x in range(1, q)}
chi = lambda a: 0 if a % q == 0 else (1 if a % q in QR else -1)
n = 12; S = np.zeros((n, n), int); S[0, 1:] = 1; S[1:, 0] = -1
for i in range(q):
    for j in range(q): S[1+i, 1+j] = chi(j - i)
Hm = np.eye(n, dtype=int) + S
Pp = np.array([[1 if Hm[p, x] != Hm[p, 0] else 0 for p in range(1, n)] for x in range(1, n)])
A = model(Pp, 1); print("Paley12 L 1 |A|", len(A), "H", exact_H(A, 4))
# Syl16 L=1 (d=15): enumerate sparse vectors only (ball: <=3 nonzeros incl base); direct filter
P16 = syl(4)
A = []
for supp in range(0, 4):
    for X in itertools.combinations(range(15), supp):
        for sg in itertools.product((-1, 1), repeat=supp):
            k = [0]*15
            for x, s_ in zip(X, sg): k[x] = s_
            if np.all(np.abs(np.array(k) @ P16) <= 1): A.append(tuple(k))
print("Syl16 L 1 |A| (supp<=3 enumeration)", len(A), "H", exact_H(A, 4))
