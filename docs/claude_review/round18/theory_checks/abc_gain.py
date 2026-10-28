"""abc-gain of unit pencils f = s*A/C on Mbar_{0,8} (anchor = infinity, integral chart x_1..x_7).
A, C: coprime products of differences (x_i - x_j) of equal degree D; B = C - s*A.
gain(f) = L.beta - F_red.beta, with (Kapranov chart, H.beta = 31, E_T.beta = 1):
  L.beta  = 31 D - sum_{|T|>=3} mu_T,  mu_T = min(nu_T(A), nu_T(C))
  F_red.beta = #{visible T : max(nu_T(A),nu_T(B),nu_T(C)) > mu_T} + sum_{interior Q-irred factors Y of B} Y.beta,
  Y.beta = 31 deg Y - sum_{|T|>=3} nu_T(Y).
Positive gain would make abc imply (L_8) on pairwise-coprime profiles (theory.md §2.5).
usage: python3 abc_gain.py Dmax
"""
import itertools, random, sys
from fractions import Fraction
import flint

N = 7
ctx = flint.fmpz_mpoly_ctx.get(('x', N), 'lex')
X = ctx.gens()
PAIRS = list(itertools.combinations(range(N), 2))
VIS = [frozenset(T) for r in range(2, N) for T in itertools.combinations(range(N), r)]
assert len(VIS) == 119
BIG = [T for T in VIS if len(T) >= 3]
rng = random.Random(12345)

def nu_mono(mono, T):  # mono: dict pair->exp
    return sum(e for (i, j), e in mono.items() if i in T and j in T)

def lagrange_low(vals, ts):
    # exact interpolation, return lowest nonzero coefficient index
    n = len(ts)
    coeffs = [Fraction(0)] * n
    for k in range(n):
        # basis poly prod_{m!=k} (t - t_m)/(t_k - t_m)
        basis = [Fraction(1)]
        den = Fraction(1)
        for m in range(n):
            if m == k: continue
            basis = [Fraction(0)] + basis
            for q in range(len(basis) - 1):
                basis[q] -= ts[m] * basis[q + 1]
            den *= (ts[k] - ts[m])
        for q in range(n):
            coeffs[q] += vals[k] * basis[q] / den
    for q, c in enumerate(coeffs):
        if c != 0: return q
    return None

CT = flint.fmpz_mpoly_ctx.get(('t',), 'lex')
TT, = CT.gens()

def nu_poly(F, deg, T):
    """order of F along L_T: lowest t-power of F(x_i = a + b_i t (i in T), x_j = r_j), two random trials."""
    best = None
    for trial in range(2):
        a = rng.randint(-10**6, 10**6)
        args = []
        for i in range(N):
            if i in T:
                args.append(CT.from_dict({(0,): a}) + rng.randint(1, 10**6) * TT)
            else:
                args.append(CT.from_dict({(0,): rng.randint(-10**9, 10**9)}))
        G = F.compose(*args, ctx=CT)
        lo = min(k[0] for k in G.to_dict().keys())
        best = lo if best is None else min(best, lo)
    return best

def nu_poly_old(F, deg, T):
    T = sorted(T)
    best = None
    for trial in range(2):
        a = rng.randint(-10**6, 10**6)
        b = {i: rng.randint(1, 10**6) for i in T}
        r = {j: rng.randint(-10**9, 10**9) for j in range(N) if j not in T}
        ts = list(range(1, deg + 2))
        vals = []
        for t in ts:
            pt = [a + b[i] * t if i in b else r[i] for i in range(N)]
            vals.append(Fraction(int(F(*pt))))
        lo = lagrange_low(vals, ts)
        best = lo if best is None else min(best, lo)
    return best

def mono_poly(mono):
    P = ctx.from_dict({(0,)*N: 1})
    for (i, j), e in mono.items():
        P = P * (X[i] - X[j]) ** e
    return P

def gain(A, C, s):
    D = sum(A.values())
    PA, PC = mono_poly(A), mono_poly(C)
    B = PC - s * PA
    if B == 0: return None
    c, facs = B.factor()
    nuB = {T: 0 for T in VIS}
    interior = 0
    for (Fk, mult) in facs:
        # pair factor?
        dk = Fk.total_degree()
        is_pair = False
        if dk == 1:
            for (i, j) in PAIRS:
                if Fk == X[i] - X[j] or Fk == X[j] - X[i]:
                    is_pair = True
                    for T in VIS:
                        if i in T and j in T: nuB[T] += mult
                    break
        if not is_pair:
            nus = {T: nu_poly(Fk, dk, T) for T in BIG}
            for T in BIG: nuB[T] += mult * nus[T]
            interior += 31 * dk - sum(nus.values())
    Lb = 31 * D
    supp = 0
    for T in VIS:
        a, cc, bb = nu_mono(A, T), nu_mono(C, T), nuB[T]
        mu = min(a, cc) if len(T) >= 3 else 0
        if len(T) >= 3: Lb -= mu
        if max(a, bb, cc) > mu: supp += 1
    return Lb - supp - interior, Lb, supp, interior

def monos(D):
    for combo in itertools.combinations_with_replacement(range(len(PAIRS)), D):
        m = {}
        for k in combo: m[PAIRS[k]] = m.get(PAIRS[k], 0) + 1
        yield m

def canon(m):
    # canonical form under S_7 (brute force)
    best = None
    for p in itertools.permutations(range(N)):
        key = tuple(sorted((min(p[i], p[j]), max(p[i], p[j]), e) for (i, j), e in m.items()))
        if best is None or key < best: best = key
    return best

if __name__ == '__main__':
    # sanity: Ptolemy and anchored Ptolemy
    print('anchored Ptolemy', gain({(0,1):1}, {(0,2):1}, 1))
    print('Ptolemy', gain({(0,1):1,(2,3):1}, {(0,2):1,(1,3):1}, 1))
    Dmax = int(sys.argv[1])
    for D in range(1, Dmax + 1):
        allm = list(monos(D))
        reps = {}
        for m in allm:
            reps.setdefault(canon(m), m)
        best = None; hist = {}
        n = 0
        for key, A in reps.items():
            for C in allm:
                if set(A) & set(C): continue
                for s in (1, -1):
                    g = gain(A, C, s)
                    if g is None: continue
                    n += 1
                    hist[g[0]] = hist.get(g[0], 0) + 1
                    if best is None or g[0] > best[0]:
                        best = (g[0], g, A, C, s)
        print(f"D={D}: A-orbits {len(reps)}, pencils {n}, max gain {best[0]} {best[1:]}")
        print("   top gains:", sorted(hist.items(), reverse=True)[:6])
        sys.stdout.flush()
