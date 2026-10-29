"""How do J9a / J9a' contain M8'?  Compare slices with the M8' box (exact set comparison)."""
import sys
import numpy as np
sys.path.insert(0, '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round28/verify_checks')
from vlib import *
n8, p8, L8, M = spec_box('syl:3:3')          # points 0..7, base 7
Mset = set(tuple(k.tolist()) for k in M)
sigM = M.astype(int).sum(axis=1)
print('M8\' |A|', len(M), ' with |k_base|<=1:', int((np.abs(sigM) <= 1).sum()))
# J9a: points 0..8, base 8; coords k_0..k_7. Slice k_8 = 0  <=>  sum_{x<8} k_x = 0.
_, _, _, A = spec_box('rsyl:4:0,1,2,3,4,5,6,7,8:3')
sl = A[A.astype(int).sum(axis=1) == 0]
# map to M8' coordinates: points 0..7 with base 7: k' = (k_0..k_6), k_7 = -sum -> automatically consistent
img = set(tuple(k[:7].tolist()) for k in sl)
inM = all(t in Mset for t in img)
sub = set(t for t in Mset if abs(sum(t)) <= 1)
print('J9a slice k_8=0:', len(sl), 'points; image in M8\':', inM, '; equals M8\' cap {|k_7|<=1}:', img == sub)
# J9a': J = (1,...,8,0), base = point 0; coords = points 1..8 (index 0..7 -> point x+1). Slice k(point 8) = 0: coord 7.
_, _, _, B = spec_box('rsyl:4:1,2,3,4,5,6,7,8,0:3')
sl2 = B[B[:, 7] == 0]
# M8' with points ordered (1,...,7, base 0)?  Build M8'' = restriction of H_16 to (1,..,7,0) (= Sylvester-8, base 0)
_, _, _, M2 = spec_box('rsyl:4:1,2,3,4,5,6,7,0:3')
M2set = set(tuple(k.tolist()) for k in M2)
img2 = set(tuple(k[:7].tolist()) for k in sl2)
print('M8 config with base 0: |A| =', len(M2), "(M8' has", len(M), ")")
print("J9a' slice k_8=0:", len(sl2), 'points; equals the M8\' box (base 0):', img2 == M2set)
