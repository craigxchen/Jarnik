"""R19-J1b (judge, unrefereed): n'(f) = #{T : f|_{D_T} constant} and Cor 2.5 bound Lambda_k - n'(f) on the 140 records.
f|_{D_T} is constant iff e(T) != 0 (value 0 or inf) or e(T) = 0 and the initial forms in_T(A), in_T(C) are proportional.
Initial form = lowest t-coefficient after x_i -> x_i0 + t x_i (i in T); exact, python-flint."""
import json, collections, sys
sys.path.insert(0, '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round19/judge_checks')
from j1_bound_vs_records import ctx, x, t, VIS, poly, inside, kf, Lam, lam_pq, OUT

def initial(p, T):
    i0 = T[0]
    sub = [x[i] if i not in T or i == i0 else x[i0] + t * x[i] for i in range(7)] + [t]
    q = p.compose(*sub)
    d = q.to_dict(); m = min(k[7] for k in d)
    return {k[:7]: v for k, v in d.items() if k[7] == m}

def proportional(a, b):
    if set(a) != set(b): return None
    k0 = next(iter(a)); r = (b[k0], a[k0])
    return r if all(b[k] * r[1] == a[k] * r[0] for k in a) else None

tot = collections.Counter()
for D in (4, 5, 6):
    for l in open(f'{OUT}/D{D}.jsonl'):
        for g, A, C, lam in [tt for tt in json.loads(l)['top'] if isinstance(tt[0], int)]:
            PA, PC = poly(A), poly(C); k = kf(A, C); npr = 0
            for T in VIS:
                if inside(A, T) != inside(C, T): npr += 1; continue
                if proportional(initial(PA, T), initial(PC, T)) is not None: npr += 1
            tot[(D, k, npr, g, Lam(k) - npr, g <= Lam(k) - npr, g == Lam(k) - npr)] += 1
out = [f"D={D} k={k} n'={n} gain={g} Cor2.5 bound={b} ok={ok} equality={eq} count={c}" for (D, k, n, g, b, ok, eq), c in sorted(tot.items())]
out.append(f"ALL_OK={all(k[5] for k in tot)} ALL_EQUALITY={all(k[6] for k in tot)} total={sum(tot.values())}")
open('/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round19/judge_checks/j1b_out.txt', 'w').write('\n'.join(out) + '\n')
print('\n'.join(out))
