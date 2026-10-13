# Lemma 5.3 count for the Kapranov foliation by lines through a point p of P^4 (exact, Z[t]).
# Kapranov model w.r.t. label 6: labels 0..4 -> e_0..e_4, label 5 -> (1,...,1), label 6 -> q, where
# q_i = (x_5 - x_i)/(x_6 - x_i). Cross-ratio of labels {a,b,c,d}: project the rational normal curve
# from the plane spanned by the other three labels; brackets [xy]_k = det[y_e, y_f, y_g, y_x, y_y].
# Leaf through q: the line q + u p. Expected: lambda(beta) = 15 - 41 = -26 for general p (local rule),
# and lambda = lambda_{f_i}(beta) = -9 for p = e_i (forgetful fibres; Lemma 4.1).
import sys, random, itertools
import flint
import lamA_exact as E
from lamA_exact import Z, ZERO, X, lam_from_br_k
def det5(M):
    # Laplace expansion along the first row (entries fmpz_poly)
    n = len(M)
    if n == 1: return M[0][0]
    s = ZERO
    for j in range(n):
        if M[0][j].is_zero(): continue
        minor = [row[:j] + row[j + 1:] for row in M[1:]]
        term = M[0][j] * det5(minor)
        s = s + term if j % 2 == 0 else s - term
    return s
def brackets(mob, p):
    N = [Z(m[1]) + Z(m[0]) * X for m in mob]; D = [Z(m[3]) + Z(m[2]) * X for m in mob]
    br = lambda a, b: N[a] * D[b] - N[b] * D[a]
    # q_i proportional to [5 i] / [6 i]  ->  q_i = [5 i] * prod_{j != i} [6 j]
    q = []
    for i in range(5):
        v = br(5, i)
        for j in range(5):
            if j != i: v = v * br(6, j)
        q.append(v)
    q = E.content_div(q)
    pts0 = [[Z(1) if j == i else ZERO for j in range(5)] for i in range(5)] + [[Z(1)] * 5, q]
    pts1 = [[ZERO] * 5 for _ in range(6)] + [[Z(c) for c in p]]          # u-derivative: only label 6 moves
    out0 = {}; out1 = {}
    for k in itertools.combinations(range(7), 4):
        rest = [a for a in range(7) if a not in k]
        base = [pts0[e] for e in rest]
        for (x, y) in itertools.permutations(k, 2):
            M0 = base + [pts0[x], pts0[y]]
            d0 = det5([list(r) for r in M0])
            # derivative: multilinear in rows; only row of label 6 depends on u
            d1 = ZERO
            rows = rest + [x, y]
            for r_i, lab in enumerate(rows):
                if lab == 6:
                    M1 = [pts0[l] if j != r_i else pts1[6] for j, l in enumerate(rows)]
                    d1 = d1 + det5([list(r) for r in M1])
            out0[(k, x, y)] = d0; out1[(k, x, y)] = d1
    return out0, out1
if __name__ == '__main__':
    random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 5)
    for name, p in [('general p', [3, -7, 2, 11, 5]), ('p = e_0 (f_O expected)', [1, 0, 0, 0, 0]),
                    ('p = e_3 (f_b1 expected)', [0, 0, 0, 1, 0])]:
        lm = lam_from_br_k(*brackets(E.int_mob(()), p))
        tot = 0; diffs = {}
        for T in itertools.combinations(range(7), 3):
            v = lam_from_br_k(*brackets(E.int_mob(T), p)); tot += v
            if v != lm: diffs[T] = v - lm
        print(name, ': lambda(mu) =', lm, ' lambda(beta) =', tot - 32 * lm, ' nonzero diffs:', len(diffs),
              sorted(set(diffs.values())), flush=True)
