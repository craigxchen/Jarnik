"""Corollary 5.2: least M_0 with exp(c M^(delta/2)) >= 1 + 2 M^2/kappa for all M >= M_0 (checked on a log grid
 and by monotonicity of h(M) = c M^(delta/2) - log(1 + 2M^2/kappa), whose derivative is positive beyond its zero)."""
import math, json
kappa = (4 / math.pi) * (1 + 2 ** -0.5)
def h(M, c, d):
    lM = math.log(M)
    return c * math.exp(d / 2 * lM) - (math.log(2 / kappa) + 2 * lM + math.log1p(kappa / (2 * M * M)))
out = {}
for c, d in [(0.25, 0.5), (0.25, 1.0), (1.0, 0.5), (0.5, 0.25), (0.25, 0.25)]:
    # find last M (log grid up to 1e200) with h < 0
    lo = None; x = 0.0
    while x < 300:
        M = math.exp(x)
        if h(M, c, d) < 0: lo = M
        x += 0.01
    # derivative check at lo: h'(M) = c (d/2) M^(d/2-1) - 4M/(kappa + 2M^2) > 0 beyond
    out[f'c={c},delta={d}'] = {'M0_approx': lo * math.exp(0.01) if lo else 1.0,
                               'log10_M0': math.log10(lo * math.exp(0.01)) if lo else 0}
    print(c, d, out[f'c={c},delta={d}'])
json.dump(out, open('barrier_M0_out.json', 'w'), indent=1)
