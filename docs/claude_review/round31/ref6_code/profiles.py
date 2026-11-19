"""(a),(b): independent construction of the record's fake profiles (binary, e=1), anchored parity rank,
number of DISTINCT nonzero anchored patterns, and the weighted Markov packing certificate
   exists j:  2^rank * (n - 8j) * C(n,j) > n * 2^n      (contradiction with (VI.1) for ANY weights).
Profiles:
  DC(q)    : Paley skew double core, round8/lw.md 5.1, full M=q+1 rows (border included)
  DCo(q)   : same, ordinary q rows only (Astra's Part III.2 / VI.1 version)
  FF(q,b)  : construct.md 2: b unflipped copies of each label of Paley-I core + one flipped copy per row
  P1(q)    : single Paley core (no flips), q nonconstant columns
  SYL(t)   : Sylvester/Walsh core, 2^t rows, 2^t-1 columns
"""
from math import comb, log2
import sys

def legendre(a, q):
    a %= q
    if a == 0:
        return 0
    return 1 if pow(a, (q - 1) // 2, q) == 1 else -1

def rank2(rows):
    piv = {}
    for row in rows:
        while row:
            l = row.bit_length()
            if l not in piv:
                piv[l] = row; break
            row ^= piv[l]
    return len(piv)

def analyse(name, S):
    """S: list of rows, each a list of +-1 (columns = primes, e=1). bits = (1+s)/2."""
    M = len(S)
    r = len(S[0])
    bits = [[(1 + s) // 2 for s in row] for row in S]
    unanch = rank2([sum(b << j for j, b in enumerate(row)) for row in bits])
    pats = {}
    for j in range(r):
        pat = tuple(bits[i][j] ^ bits[0][j] for i in range(1, M))
        if any(pat):
            pats[pat] = pats.get(pat, 0) + 1
    P = list(pats)
    n = len(P)
    rows = [sum(P[k][i] << k for k in range(n)) for i in range(M - 1)]
    rk = rank2(rows)
    # certificate
    cert = None
    for j in range(0, n // 8 + 1):
        if n - 8 * j <= 0:
            break
        if (1 << rk) * (n - 8 * j) * comb(n, j) > n * (1 << n):
            cert = j
            break
    print(f"{name:12s} M={M:4d} cols={r:5d} unanch.rank={unanch:4d} distinct anchored patterns n={n:4d} "
          f"anchored rank={rk:4d} rank/n={rk/n:.3f}  certificate j={cert}")
    return rk, n, cert

CHECK=True
def paley_K(q):
    # rows/cols indexed 0 (border), 1..q  <-> F_q element x-1
    M = q + 1
    K = [[0] * M for _ in range(M)]
    for x in range(M):
        for y in range(M):
            if x == 0 and y == 0:
                S = 0
            elif x == 0:
                S = -1
            elif y == 0:
                S = 1
            else:
                S = legendre((y - 1) - (x - 1), q)
            K[x][y] = (1 if x == y else 0) + S
    # checks: skew-Hadamard
    if CHECK and q < 120:
        for x in range(M):
            for y in range(M):
                ip = sum(K[x][k] * K[y][k] for k in range(M))
                assert ip == (M if x == y else 0)
    return K

def double_core(q, border=True):
    K = paley_K(q)
    M = q + 1
    cols = []
    for a in range(1, M):
        c1 = [K[x][a] for x in range(M)]
        c2 = list(c1)
        c2[0] = -c2[0]
        c2[a] = -c2[a]
        cols += [c1, c2]
    S = [[c[x] for c in cols] for x in range(M)]
    return S if border else S[1:]

def paley_H(q):
    # construct.md Paley-I: H(inf,.)=1, H(i,inf)=-1, H(i,j)=delta_ij+chi(j-i); then rows i in F_q negated
    M = q + 1
    H = [[0]*M for _ in range(M)]
    for x in range(M):
        for y in range(M):
            if x == 0: H[x][y] = 1
            elif y == 0: H[x][y] = -1
            else: H[x][y] = (1 if x == y else 0) + legendre((y-1)-(x-1), q)
    Hn = [[H[x][y] * (1 if x == 0 else -1) for y in range(M)] for x in range(M)]
    for x in range(M):
        assert Hn[x][0] == 1
    if q < 120:
        for x in range(M):
            for y in range(M):
                assert sum(Hn[x][k]*Hn[y][k] for k in range(M)) == (M if x == y else 0)
    return Hn

def flipped(q, b):
    Hn = paley_H(q)
    M = q + 1
    cols = []
    for a in range(1, M):
        for _ in range(b):
            cols.append([Hn[x][a] for x in range(M)])
    astar = 1
    for x in range(M):
        a = x if x != 0 else astar
        c = [Hn[y][a] for y in range(M)]
        c[x] = -c[x]
        cols.append(c)
    return [[c[x] for c in cols] for x in range(M)]

def single_paley(q):
    Hn = paley_H(q)
    M = q + 1
    return [[Hn[x][a] for a in range(1, M)] for x in range(M)]

def sylvester(t):
    M = 2 ** t
    return [[(-1) ** bin(x & a).count('1') for a in range(1, M)] for x in range(M)]

if __name__ == "__main__":
    for q in (167, 191, 251):
        analyse(f"DCo({q})", double_core(q, border=False))
        analyse(f"DC({q})", double_core(q, border=True))
    primes=[q for q in range(7, 420, 4) if all(q % p for p in range(2, int(q ** 0.5) + 1))]
    for fam, fn in (("DCo", lambda q: double_core(q, False)), ("DC", lambda q: double_core(q, True)),
                    ("FF3", lambda q: flipped(q, 3)), ("P1", single_paley)):
        first=None
        for q in primes:
            rk,n,cert=analyse(f"{fam}({q})", fn(q)) if q in (7,11,19,23,43) else (None,None,None)
        # threshold scan quietly
        import io, contextlib
        for q in primes:
            buf=io.StringIO()
            with contextlib.redirect_stdout(buf):
                rk,n,cert=analyse(f"{fam}({q})", fn(q))
            if cert is not None:
                print("  first certificate:", buf.getvalue().strip()); break
        else:
            print("  no certificate up to q=", primes[-1])
        print()
    for t in range(3, 9):
        analyse(f"SYL({t})", sylvester(t))
