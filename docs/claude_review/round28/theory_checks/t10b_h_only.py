"""T10b: h = max sigma on the Sylvester T=1 box (IP, long time limit), to decide rigidity where t10 timed out.
Usage: t10b_h_only.py k L timelimit"""
import sys
exec(open('/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round28/theory_checks/t10_bco_ip.py').read().split("t0 = time.time()")[0])
import time
t0 = time.time()
r = solve_h()
print("Sylvester n=%d L=%d: h incumbent %s, dual bound %s, status %s (%.0fs)" % (n, L, -r.fun if r.x is not None else None,
      getattr(r, 'mip_dual_bound', None), r.status, time.time() - t0), flush=True)
