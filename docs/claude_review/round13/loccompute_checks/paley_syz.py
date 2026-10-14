import itertools, numpy as np, time
from loclib import Ring
from cyclib import CycRing, ideal_dims_cyc
from cyc_run import dir_orbit_reps, image_sizes, powgen
import models
d, cols, Ls = models.paley12(2); A = models.box(d, cols, Ls)
ring = Ring(d, 10); cr = CycRing(ring, 11)
V = dir_orbit_reps(d, 1); sz = image_sizes(A, V, -22, 22)
q5 = [powgen(ring, v, 5) for v, s in zip(V, sz) if s == 5]
q5 = [g for g in q5 if g.any()]
print("quintic width gens (orbit reps):", len(q5), [ [v for v,s in zip(V,sz) if s==5] ])
dims = ideal_dims_cyc(cr, {5: np.array(q5)}, 10, verbose=False)
print("dims of ideal (12 quintics) in R, deg 0..10:", dims)
print("N_i:", ring.N[:11])
print("free-module count 12*N_{i-5}:", [12*ring.N[i-5] if i>=5 else 0 for i in range(11)])
