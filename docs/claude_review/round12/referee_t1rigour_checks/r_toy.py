# Referee: Theorem 6.11(a) toy model, exact. A = {k in {-1,0,1}^6 : |k1+k2+k3+k4| <= 1}.
import itertools, flint
A = [k for k in itertools.product((-1,0,1), repeat=6) if abs(sum(k[:4])) <= 1]
print("|A| =", len(A))
def monos(s, e):
    return [m for deg in range(e+1) for m in itertools.product(range(deg+1), repeat=s) if sum(m) == deg]
def ev(k, m):
    v = 1
    for a, b in zip(k, m): v *= a**b
    return v
def rank_rows(rows):
    return flint.fmpz_mat(rows).rank()
cache = {}
for m in range(6, 10):
    ms = monos(6, m-1)
    Mt = [[ev(k, mo) for k in A] for mo in ms]   # rows = monomials (functions on A)
    r = rank_rows(Mt)
    dimV = len(A) - r
    res = []
    for u in ((1,3,7,19,0,0), (1,3,7,19,5,0), (1,3,7,19,43,101)):
        row = [sum(a*b for a, b in zip(k, u))**m for k in A]
        r2 = rank_rows(Mt + [row])
        res.append((u, r2 == r))   # True <=> <k,u>^m in Poly_{<m}|_A <=> Phi_c(u)=0 for all c in V(A,m)
    print(f"m={m}: dim V(A,m) = {dimV}; common zero of all Phi_c:", res)
# zonotope inside A with sum L_i = 8 (generators e1-e2, e3-e4, e5, e6, each L=2)
k0 = (-1,1,-1,1,-1,-1)
gens = [(1,-1,0,0,0,0),(0,0,1,-1,0,0),(0,0,0,0,1,0),(0,0,0,0,0,1)]
Aset = set(A)
inside = all(tuple(k0[i] + sum(t[g]*gens[g][i] for g in range(4)) for i in range(6)) in Aset
             for t in itertools.product(range(3), repeat=4))
print("zonotope k0 + sum_{t in [0,2]^4} t_i g_i inside A:", inside, " => zeta(A) >= 8, so V(A,m) has a character certificate for every m <= 8")
