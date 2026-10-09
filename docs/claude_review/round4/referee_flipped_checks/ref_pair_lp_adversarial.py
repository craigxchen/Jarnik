"""Adversarial (hill-climbing) search over assignments x -> a(x) with capacity <= B for the largest
pair-LP optimum rho = max W_F/W (see ref_pair_lp_sweep.py).  Exact rational LP; certificate verified."""
import random, sys
sys.path.insert(0, '.')
from ref_pair_lp import paley, sylv
from ref_pair_lp_sweep import solve

def climb(H, B, rng, iters):
    M = len(H)
    # random start respecting capacity
    cnt = {}; assign = []
    for x in range(M):
        ch = [a for a in range(1, M) if cnt.get(a, 0) < B]; a = rng.choice(ch); assign.append(a); cnt[a] = cnt.get(a, 0) + 1
    best, _ = solve(H, assign)
    for it in range(iters):
        x = rng.randrange(M); a = rng.randrange(1, M)
        if a == assign[x] or cnt.get(a, 0) >= B: continue
        new = list(assign); new[x] = a
        v, _ = solve(H, new)
        if v >= best:
            cnt[assign[x]] -= 1; cnt[a] = cnt.get(a, 0) + 1; assign = new; best = v
    return best, assign

if __name__ == '__main__':
    rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 3)
    for name, H, B, starts, iters in [('Walsh-8', sylv(3), 2, 6, 60), ('Walsh-8', sylv(3), 3, 6, 60),
                                       ('Paley-12', paley(11), 2, 4, 80), ('Paley-12', paley(11), 3, 4, 80),
                                       ('Paley-12', paley(11), 5, 3, 80), ('Paley-12', paley(11), 4, 3, 80)]:
        res = [climb(H, B, rng, iters) for _ in range(starts)]
        v, a = max(res, key=lambda t: t[0])
        print(f"{name} B={B}: best rho found {v} = {float(v):.4f}  assignment {a}", flush=True)
