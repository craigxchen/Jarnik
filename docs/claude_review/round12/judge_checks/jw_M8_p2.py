"""Second-prime recheck of dim (J_w)_i for M8 (the degree-6 gap 616 < 637 is a mod-p statement)."""
import sys; sys.path.insert(0, '.')
import jw_test
jw_test.PR = 1000000007
import jproj_test; jproj_test.PR = 1000000007
from jw_test import widths, sylvester8, eval_rank_profile
from jproj_test import powl, jdims
import math
A = sylvester8(2); d = 7; m = 8
cache = {}; g = {}
for v, h in widths(A, d, 2, m-1):
    g.setdefault(h+1, []).append(powl(v, h+1, d, cache))
r = eval_rank_profile(A, d, m)
H = [r[0]] + [r[i]-r[i-1] for i in range(1, m+1)]
s = f"p=1000000007: dim J_w = {jdims(g, d, m)}; dim J(A) = {[math.comb(i+d-1,d-1)-H[i] for i in range(m+1)]}"
print(s); open('jw_M8_p2_out.txt','w').write(s+"\n")
