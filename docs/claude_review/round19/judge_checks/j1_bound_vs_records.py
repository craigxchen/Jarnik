"""R19-J1 (judge, unrefereed): test theory.md Thm 3.3 on every extremal record of search.md.

For each record [gain, A, C, lambda] with gain > -40 in search_checks/out/D{4,5,6}.jsonl:
  * e = edges(A) - edges(C) (chart indices 0..6, anchor 8 at infinity);
  * k(f) = |I| + [some vertex sum != 0]  (referee_theory m1);
  * N(u) = #{T visible : a_T != c_T};
  * n_1 = #{T visible : nu_T(B) > mu_T},  B = q C - p A, nu_T exact by substitution x_i -> x_i0 + t x_i (i in T);
  * Thm 3.3 bound  Lambda_k - N - n_1, Lambda_4 = 16, Lambda_k = 32 - 2^(8-k).
All arithmetic exact (python-flint fmpz_mpoly). Writes j1_out.txt.
"""
import json, itertools, collections, sys
from fractions import Fraction
import flint

OUT = '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round19/search_checks/out'
names = tuple(f'x{i}' for i in range(7)) + ('t',)
ctx = flint.fmpz_mpoly_ctx.get(names, 'lex')
X = ctx.gens()
x, t = X[:7], X[7]
VIS = [T for r in range(2, 7) for T in itertools.combinations(range(7), r)]
assert len(VIS) == 119


def poly(edges):
    p = ctx.from_dict({(0,) * 8: 1})
    for i, j in edges:
        p = p * (x[i] - x[j])
    return p


def nu(p, T):
    i0 = T[0]
    sub = [x[i] if i not in T or i == i0 else x[i0] + t * x[i] for i in range(7)] + [t]
    q = p.compose(*sub)
    return min(m[7] for m in q.monoms())


def inside(edges, T):
    s = set(T)
    return sum(1 for i, j in edges if i in s and j in s)


def lam_pq(l):
    return (l, 1) if isinstance(l, int) else (l[0], l[1])


def kf(A, C):
    e = collections.Counter()
    for i, j in A:
        e[(min(i, j), max(i, j))] += 1
    for i, j in C:
        e[(min(i, j), max(i, j))] -= 1
    e = {k: v for k, v in e.items() if v}
    I = {i for k in e for i in k}
    vs = collections.Counter()
    for (i, j), v in e.items():
        vs[i] += v; vs[j] += v
    return len(I) + (1 if any(vs.values()) else 0)


def Lam(k):
    return 16 if k == 4 else 32 - 2 ** (8 - k)


def main():
    lines = []
    tot = collections.Counter()
    for D in (4, 5, 6):
        recs = []
        for l in open(f'{OUT}/D{D}.jsonl'):
            r = json.loads(l)
            recs += [tt for tt in r['top'] if isinstance(tt[0], int)]
        for g, A, C, l in recs:
            p, q = lam_pq(l)
            B = q * poly(C) - p * poly(A)
            k = kf(A, C)
            N = sum(1 for T in VIS if inside(A, T) != inside(C, T))
            n1 = sum(1 for T in VIS if nu(B, T) > min(inside(A, T), inside(C, T)))
            bound = Lam(k) - N - n1
            ok = g <= bound
            tot[(D, k, N, n1, g, bound, ok)] += 1
        lines.append(f'D={D}: {len(recs)} records')
    for key, v in sorted(tot.items()):
        D, k, N, n1, g, bound, ok = key
        lines.append(f'  D={D} k(f)={k} N={N} n_1={n1} gain={g} Thm3.3 bound={bound} ok={ok} equality={g == bound} count={v}')
    allok = all(key[-1] for key in tot)
    alleq = all(key[4] == key[5] for key in tot)
    lines.append(f'ALL_BOUND_OK={allok}  ALL_EQUALITY={alleq}  total={sum(tot.values())}')
    open('/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round19/judge_checks/j1_out.txt', 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == "__main__":
    main()
