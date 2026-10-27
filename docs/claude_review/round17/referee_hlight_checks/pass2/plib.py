# pass2 referee library: sdc_q built from fake.md section 1 description:
# rows F_q u {inf}; S_a = {x in F_q u {inf}: K(x,a) = -1} for a Paley-I skew core, S'_a = S_a xor {inf, a}.
# We take K(x,y) = chi(y-x) for x != y in F_q, K(x,x) = 1, K(inf,y) = 1, K(x,inf) = -1 (skew-Hadamard, Paley I);
# the opposite convention chi(x-y) is the reflection x -> -x (checked: both give identical sigma).
def chi_table(q):
    t = [0]*q
    for x in range(1, q): t[(x*x) % q] = 1
    return [0] + [1 if t[u] else -1 for u in range(1, q)]
def patterns(q, conv=+1):
    ch = chi_table(q); INF = q
    def K0(x, y):
        if x == INF: return 1
        if y == INF: return -1
        if x == y: return 1
        return ch[(conv*(y - x)) % q]
    # exact.md Thm 4.1 normalisation: K -> D K D with D = diag(K(x, base)), base = INF, so that K(x, base) = 1
    def K(x, y):
        return K0(x, INF) * K0(x, y) * K0(y, INF)
    pats = []
    for a in range(q):
        S = frozenset(x for x in list(range(q)) + [INF] if K(x, a) == -1)
        Sp = S.symmetric_difference({INF, a})
        pats += [S, Sp]
    return pats
def N_def(lam, pats):
    # lam: dict/list indexed 0..q (q = inf); returns sum_S |lam(S)|
    return sum(abs(sum(lam[x] for x in S)) for S in pats)
def N_formula(lam, q):
    ch = chi_table(q); f = lam[:q]; F = sum(f)
    tot = 0
    for a in range(q):
        c = sum(f[t]*ch[(a - t) % q] for t in range(q))
        tot += max(abs(F + f[a]), abs(c))
    return tot
