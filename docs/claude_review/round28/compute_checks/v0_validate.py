"""v0: validate slib against the round-27 killer's independent library (klib: exact over Q, both signs of every shift)
on small random slab boxes, all degrees.  Also checks M_e(A_{-l}) = M_e(A_l) implicitly (slib uses one sign)."""
import sys, random, time
sys.path.insert(0, '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round28/compute_checks')
sys.path.insert(1, '/home/user/Jarnik/docs/claude_review/round27/kill_interp_checks')
import slib
import klib
out = []
def pr(s):
    print(s, flush=True); out.append(s)
rng = random.Random(2028)
cases = [(4, []), (3, [((0, 1), 1)])]   # cube d=3, hexagon
for t in range(10):
    n = rng.choice([4, 5])
    r = rng.randint(1, n + 2)
    pats = set()
    while len(pats) < r:
        S = tuple(x for x in range(n) if rng.random() < 0.5)
        if 0 < len(S) < n:
            pats.add(slib.canon(n, S))
    L = rng.choice([1, 2])
    cases.append((n, [(S, L) for S in sorted(pats)]))
bad = 0
for n, slabs in cases:
    d = n - 1
    A = slib.build_box(n, slabs)
    Ak = klib.build_box(n, slabs)
    assert [tuple(int(v) for v in k) for k in A] == [tuple(k) for k in Ak], "box mismatch"
    for p in (slib.P1, slib.P2):
        H, _ = slib.hilbert(A, d, p)
        rho = len(H) - 1
        degs = list(range(rho + 1))
        acc = slib.run_shifts(A, d, H, slib.Lpair(d), p, degs)
        mine = [(acc.rank[e] if not acc.full[e] else H[e], H[e]) for e in degs]
        if p == slib.P1:
            ex = [klib.shift_dims_Q(Ak, klib.Lpair(d), e) for e in degs]
            exs = [(b, a) for a, b in ex]
            ok = all(m == x for m, x in zip(mine, exs))
            bad += (not ok)
            pr("n=%d slabs=%s |A|=%d H=%s  slib(sum,H)=%s  klib_Q(sum,H)=%s  %s" % (
                n, slabs, len(A), H, mine, exs, "AGREE" if ok else "DISAGREE"))
        else:
            pr("   P2: H=%s sum=%s" % (H, mine))
pr("disagreements: %d" % bad)
open('/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round28/compute_checks/v0_validate_out.txt', 'w').write("\n".join(out) + "\n")
