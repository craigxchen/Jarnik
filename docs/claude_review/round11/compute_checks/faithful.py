"""For every screened instance (scanN_out.txt lines), compute zeta(A) exactly (zono.py) and the exact rho(A)
(fmpz rank), and list the faithful windows  max(zeta, 2*rmax) < m <= rho  [II]  and  max(zeta, rmax) < m <= rho [II']."""
import sys, ast, re
from zono import zeta
from lpzeta import lp_bound
import signal
from models import rho
from core import build_A
fn = sys.argv[1]
seen = set()
for line in open(fn):
    m = re.match(r"(\d+) (\(.*\)) T=(\d+) L=(\d+) \|A\|=(\d+) rmax=(\d+) rho=(\d+) (\S+)", line.strip())
    if not m:
        continue
    n = int(m.group(1)); prof = ast.literal_eval(m.group(2)); T = int(m.group(3)); L = int(m.group(4))
    rmax = int(m.group(6)); rho_p = int(m.group(7))
    A = build_A(n, T, [(S, L) for S in prof])
    key = (n, tuple(sorted(A)))
    if key in seen:
        continue
    seen.add(key)
    rh = rho(A)
    lb, _ = lp_bound(n, T, prof, L)
    if lb < rh:
        z, Z = None, None
        zstr = "<=%s(LP)" % lb
        zub = int(lb)
    else:
        def handler(sig, frm): raise TimeoutError
        signal.signal(signal.SIGALRM, handler); signal.alarm(int(sys.argv[2]) if len(sys.argv) > 2 else 600)
        try:
            z, Z = zeta(n, T, prof, L, rh); zstr = str(z); zub = z
        except TimeoutError:
            z = None; zstr = "undetermined(LP<=%s, timeout)" % lb; zub = rh
        signal.alarm(0)
    z = zub
    winII = [mm for mm in range(max(z, 2 * rmax) + 1, rh + 1)]
    winII_ = [mm for mm in range(max(z, rmax) + 1, rh + 1)]
    print(n, prof, "T=%d L=%d |A|=%d rmax=%d rho=%d zeta=%s  faithful II:%s II':%s" % (T, L, len(A), rmax, rh, zstr, winII, winII_), flush=True)
