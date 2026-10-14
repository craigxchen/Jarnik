"""Syzygies of the column width powers y_p^{2L+1} in R = S/(x^3) versus the Froberg-type expectation
[(1+t+t^2)^d (1 - t^{2L+1})^P]_+ (positive part).  Mod p (evidence)."""
import numpy as np, itertools
from locrun import *
def run(name, d, cols, Ls, D):
    ring = Ring(d, D)
    e = 2 * Ls[0] + 1
    g = {e: [width_gens(ring, [(tuple(c), e)])[e][0] for c in cols]}
    dims = ideal_dims_modp(ring, g, D)
    N = np.array(ring.N[:D + 1])
    ser = np.array([1] + [0] * D)
    for _ in cols:
        ser = np.convolve(ser, [1] + [0] * (e - 1) + [-1])[:D + 1]
    exp_ = np.convolve(N, ser)[:D + 1] if False else None
    # expected codim = coefficients of (1+t+t^2)^d * (1-t^e)^P, truncated at first nonpositive
    c = np.array([1]); 
    for _ in range(d): c = np.convolve(c, [1, 1, 1])
    for _ in cols: c = np.convolve(c, [1] + [0] * (e - 1) + [-1])
    c = c[:D + 1]; expc = []
    for v in c:
        if v <= 0 or (expc and expc[-1] == 0): expc.append(0)
        else: expc.append(int(v))
    expdim = [int(ring.N[i]) - expc[i] for i in range(D + 1)]
    print(f"{name}: e={e} P={len(cols)}\n   actual dim (y^e) in R: {dims}\n   Froberg-expected     : {expdim}", flush=True)
run("M8 (Syl(8,2))", *models.sylvester(3, 2), 10)
run("M8' (Syl(8,3))", *models.sylvester(3, 3), 13)
run("Paley-12 L=2 (11 columns only)", *models.paley12(2), 9)
