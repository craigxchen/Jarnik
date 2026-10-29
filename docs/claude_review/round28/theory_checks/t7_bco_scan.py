"""T7: base-coordinate obstruction (BCO) data on the round-27 boxes: h = max sigma, rigidity (R) of the top level w.r.t.
pair moves e_x - e_y, and the reach (from t1 / k11 outputs, or recomputed mod p here for random boxes when asked)."""
import sys, random, re
sys.path.insert(0, '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round28/theory_checks')
from tlib import *
def k11_box(seed, n):
    rng = random.Random(seed); d = n - 1
    r = rng.randint(n - 1, 2 * n); pats = set()
    while len(pats) < r:
        S = tuple(x for x in range(n) if rng.random() < 0.5)
        if 0 < len(S) < n: pats.add(canon(n, S))
    L = rng.choice([1, 2])
    return n, [(S, L) for S in sorted(pats)]
def bco(n, slabs):
    d = n - 1
    A = build_box(n, slabs); As = set(A)
    h = max(sum(k) for k in A)
    top = [k for k in A if sum(k) == h]
    nonrigid = 0
    for k in top:
        for x in range(d):
            for y in range(d):
                if x != y:
                    kk = list(k); kk[x] += 1; kk[y] -= 1
                    if tuple(kk) in As: nonrigid += 1
    return len(A), h, len(top), nonrigid
rows = []
for name in ('syl8', 'sdc7'):
    for L in (1, 2, 3):
        n, pats = syl(3) if name == 'syl8' else sdc(7)
        nA, h, nt, nr = bco(n, [(S, L) for S in pats])
        print("%s L=%d |A|=%d h=%d 2h=%d |top|=%d nonrigid pair-moves=%d -> BCO predicts reach<=2h=%d if (N)" % (name, L, nA, h, 2*h, nt, nr, 2*h) if nr == 0 else
              "%s L=%d |A|=%d h=%d 2h=%d |top|=%d nonrigid pair-moves=%d -> (R) fails" % (name, L, nA, h, 2*h, nt, nr), flush=True)
# random boxes of k11 (reach from k11 outputs, all-cube shifts; k8 = L_pair) 
k8 = {}
for fn in ['/home/user/Jarnik/docs/claude_review/round27/kill_interp_checks/k8_random_100_15_6_out.txt',
           '/home/user/Jarnik/docs/claude_review/round27/kill_interp_checks/k8_random_1_12_7_out.txt']:
    for line in open(fn):
        m = re.match(r"seed=(\d+) n=(\d+) .* rho=(\d+) .* criterion reach=(\d+) zeta(>=|=)(\d+)", line)
        if m: k8[(int(m.group(1)), int(m.group(2)))] = (int(m.group(3)), int(m.group(4)), m.group(5) + m.group(6))
for (sd, nn), (rho, reach, z) in sorted(k8.items(), key=lambda t: (t[0][1], t[0][0])):
    n, slabs = k11_box(sd, nn)
    nA, h, nt, nr = bco(n, slabs)
    print("random seed=%d n=%d L=%d |A|=%d rho=%d reach(L_pair)=%d zeta%s | h=%d 2h=%d |top|=%d nonrigid=%d %s" % (
        sd, nn, slabs[0][1], nA, rho, reach, z, h, 2*h, nt, nr,
        ("BCO-rigid: reach==2h" if reach == 2*h else "BCO-rigid: reach!=2h") if nr == 0 else "not rigid"), flush=True)
