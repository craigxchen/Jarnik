"""Is the classifying point separable?  Cross-ratios of 4 points of the cluster, and ratios z_a/z_b,
as elements of F_4(t): f = num/den is in F_4(t^2) (inseparable) iff num' den - num den' = 0."""
import pickle, itertools, sys
from f24 import *
def der4(P):
    return (deriv2(P[0]), deriv2(P[1]))
for (m, L) in [(8, 16), (16, 48), (32, 128), (32, 129)]:
    S = pickle.load(open(f'cluster_m{m}_L{L}.pkl', 'rb'))
    z, forms = S['z'], S['forms']
    nz_cr = tot_cr = 0
    for a, b, c, dd in list(itertools.combinations(forms, 4))[:60]:
        num = mul4(add4(z[a], z[c]), add4(z[b], z[dd]))
        den = mul4(add4(z[a], z[dd]), add4(z[b], z[c]))
        w = add4(mul4(der4(num), den), mul4(num, der4(den)))
        tot_cr += 1; nz_cr += (not iszero4(w))
    nz_r = tot_r = 0
    for a, b in list(itertools.combinations(forms, 2))[:60]:
        w = add4(mul4(der4(z[a]), z[b]), mul4(z[a], der4(z[b])))
        tot_r += 1; nz_r += (not iszero4(w))
    # is Y a square?  Y' = 0 iff L even
    print(f"m={m} L={L}: cross-ratios with nonzero derivative {nz_cr}/{tot_cr}; ratios z_a/z_b separable {nz_r}/{tot_r}; Y'=0: {L%2==0}")
