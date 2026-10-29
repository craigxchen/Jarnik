"""zeta for one box, classes 'pair' and 'all'.  Exact pipeline:
  B = exact rational LP-dual bound (record LP);  for t = floor(B), floor(B)-1, ...: exact DFS for a zonotope of order
  >= t over the candidates allowed by the dual reduction sum_v x_v (c_v - 1) <= B - t (c_v = sum_lam y_lam|lam(v)| >= 1);
  the first feasible t is zeta, witness re-verified exactly; 'exact' iff every DFS above it completed.
Usage: run_zeta.py SPEC [timelimit_per_t].  Appends to zeta.jsonl."""
import sys, json, time, math
sys.path.insert(0, '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round28/compute_checks')
import slib, zono
from run_box import parse, OUT

spec = sys.argv[1]
tl = int(sys.argv[2]) if len(sys.argv) > 2 else 900
n, pats, L = parse(spec)
d = n - 1
slabs = [(S, L) for S in pats]
A = slib.build_box(n, slabs)
rec = dict(spec=spec, n=n, d=d, r=len(pats), L=L, size=int(len(A)), method='LPdual+DFS')
for cls in ('pair', 'all'):
    t0 = time.time()
    V, lams, bs = zono.candidates(n, slabs, cls)
    B, f, yy = zono.lp_bound(V, lams, bs)
    complete_above = True
    found = None
    log = []
    for t in range(math.floor(B), 0, -1):
        w, complete, nodes = zono.dfs_feasible(n, slabs, V, lams, bs, yy, B, t, time_limit=tl)
        log.append((t, w is not None, complete, nodes))
        print(spec, cls, 'target', t, 'found' if w else 'none', 'complete' if complete else 'TIMEOUT', nodes,
              '%.0fs' % (time.time() - t0), flush=True)
        if w is not None:
            found = (t, w); break
        if not complete:
            complete_above = False
    if found is None:
        rec[cls] = dict(ncand=int(len(V)), lp=str(B), witness=None, log=log); continue
    t, (k0, gens) = found
    chk = zono.verify_exact(A, k0, gens, d)
    ok = chk['support_in_A'] and chk['all_zonotope_points_in_A'] and chk['in_V_t'] and chk['not_in_V_t1']
    rec[cls] = dict(ncand=int(len(V)), lp=str(B), lp_floor=int(math.floor(B)), witness_order=chk['order'],
                    witness_ok=ok, k0=k0, gens=gens, sup=chk['sup'], exact=bool(ok and complete_above),
                    upper=('LP' if chk['order'] == math.floor(B) else ('DFS-exhaustive' if complete_above else 'open')),
                    log=log, secs=round(time.time() - t0))
    print(spec, cls, rec[cls], flush=True)
with open(OUT + 'zeta.jsonl', 'a') as fh:
    fh.write(json.dumps(rec) + "\n")
