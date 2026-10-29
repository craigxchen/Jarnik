"""Driver: one slab box, one prime.  Usage:
  run_box.py SPEC PRIME [cube] [cubemax]
SPEC: 'syl:k:L', 'sdc:q:L', 'had:q:L', 'rsyl:k:J:L' / 'rsdc:q:J:L' / 'rhad:q:J:L' (J = comma list of point indices,
      last = base), 'rand:n:seed' (round-27 killer law), 'randL:n:seed:L' (same profile law, L forced).
Computes H_A (mod PRIME), reach w.r.t. L_pair (first e with dim sum_l M_e(A_l) < H_A(e); rho if none),
and with 'cube' the reach w.r.t. all {-1,0,1}^d shifts (continuing from the L_pair accumulators; only the degrees
from the L_pair failure on; cubemax = max number of extra shift classes, default all).
Appends one JSON line to results_<PRIME>.jsonl."""
import sys, json, time
sys.path.insert(0, '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round28/compute_checks')
import numpy as np
import slib, models

OUT = '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round28/compute_checks/'

def parse(spec):
    f = spec.split(':')
    if f[0] in ('syl', 'sdc', 'had'):
        M, pats = {'syl': models.syl, 'sdc': models.sdc, 'had': models.had_paley}[f[0]](int(f[1]))
        return M, pats, int(f[2])
    if f[0] in ('rsyl', 'rsdc', 'rhad'):
        M, pats = {'rsyl': models.syl, 'rsdc': models.sdc, 'rhad': models.had_paley}[f[0]](int(f[1]))
        J = [int(x) for x in f[2].split(',')]
        n, p2 = models.restrict(M, pats, J)
        return n, p2, int(f[3])
    if f[0] == 'rand':
        n = int(f[1]); r, L, pats = models.random_profile(n, int(f[2]))
        return n, pats, L
    if f[0] == 'randL':
        n = int(f[1]); r, L, pats = models.random_profile(n, int(f[2]))
        return n, pats, int(f[3])
    raise ValueError(spec)

def main():
    spec = sys.argv[1]; p = int(sys.argv[2])
    do_cube = len(sys.argv) > 3 and sys.argv[3] == 'cube'
    cubemax = int(sys.argv[4]) if len(sys.argv) > 4 else None
    n, pats, L = parse(spec)
    d = n - 1
    slabs = [(S, L) for S in pats]
    A = slib.build_box(n, slabs)
    t0 = time.time()
    log = lambda s: print(s, flush=True)
    log("%s p=%d n=%d d=%d r=%d L=%d |A|=%d" % (spec, p, n, d, len(pats), L, len(A)))
    H, _ = slib.hilbert(A, d, p)
    rho = len(H) - 1
    log("  H=%s rho=%d (%.0fs)" % (H, rho, time.time() - t0))
    degs = list(range(rho))
    acc = slib.run_shifts(A, d, H, slib.Lpair(d), p, degs, log=log, tag='pair')
    rp = slib.reach_from(acc, degs, rho)
    dims_pair = {e: (min(acc.rank[e], H[e]), H[e]) for e in degs}
    log("  reach_pair=%d dims(first failure)=%s (%.0fs)" % (rp, dims_pair.get(rp), time.time() - t0))
    rec = dict(spec=spec, p=p, n=n, d=d, r=len(pats), L=L, pats=[list(S) for S in pats], size=int(len(A)), H=H,
               rho=rho, reach_pair=rp, fail_pair=dims_pair.get(rp), t_pair=round(time.time() - t0))
    if do_cube and rp < rho:
        # continue with all cube shift classes not in L_pair, degrees rp..rho-1
        pairset = set(slib.Lpair(d))
        extra = [l for l in slib.Lcube(d) if l not in pairset]
        # degree by degree from the L_pair failure: only the lowest unfilled degree matters for the reach
        ntotal = len(extra)
        if cubemax is not None:
            sz = [(slib.shift_set(A, l).shape[0], l) for l in extra]
            sz.sort(key=lambda t: -t[0])
            extra = [l for s_, l in sz[:min(cubemax, len(extra))]]
        nextra = len(extra)
        e = rp
        last = None
        while e < rho:
            if acc.full[e]:
                e += 1; continue
            acc2 = slib.Acc(d, [e], {e: H[e]}, p)
            acc2.basis[e] = acc.basis[e]; acc2.rank[e] = acc.rank[e]
            acc2 = slib.run_shifts(A, d, H, extra, p, [e], acc=acc2, log=log, tag='cube e=%d' % e)
            last = (min(acc2.rank[e], H[e]), H[e])
            if not acc2.full[e]:
                break
            e += 1
        rc = e
        rec.update(reach_cube=rc, fail_cube=(last if rc < rho else None),
                   cube_classes=nextra + len(pairset), cube_complete=(nextra == ntotal))
        log("  reach_cube=%d (%.0fs)" % (rc, time.time() - t0))
    elif do_cube:
        rec.update(reach_cube=rp, fail_cube=None, cube_classes=0, cube_complete=True)
    rec['t_total'] = round(time.time() - t0)
    with open(OUT + 'results_%d.jsonl' % p, 'a') as fh:
        fh.write(json.dumps(rec) + "\n")

if __name__ == '__main__':
    main()
