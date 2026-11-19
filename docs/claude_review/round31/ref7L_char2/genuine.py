# (d) Exact F_4[t] instances of the genuine VII.6 model, and of the referee's dimension-4 variant.
import itertools, sys
from gf import *
W = 2  # theta
def tpow(L): return [0] * L + [1]
def elemF2(bits): return f2poly_to_f4(bits)
def build(xs_bits, L, hyperplanes):
    """xs_bits: F_2[t] elements (int bitmasks) spanning V; hyperplanes: list of lists of coefficient vectors."""
    Y = p_add(tpow(L), [W])
    marks = []
    for H in hyperplanes:
        z = [1]
        for c in H:
            v = 0
            for ci, xb in zip(c, xs_bits):
                if ci: v ^= xb
            assert v != 0
            z = p_mul(z, p_add(Y, elemF2(v)))
        marks.append(z)
    return marks
def analyse(marks, label, L, B, slope):
    norms = [p_mul(z, p_conj(z)) for z in marks]
    assert all(nm == norms[0] for nm in norms), "norms differ"
    assert all(c in (0, 1) for c in norms[0]), "norm not in F_2[t]"
    assert len({tuple(z) for z in marks}) == len(marks), "not distinct"
    D0 = deg(marks[0]); assert all(deg(z) == D0 for z in marks)
    g = marks[0]
    for z in marks[1:]: g = p_gcd(g, z)
    dg = deg(g)
    D = D0 - dg
    md = max(deg(p_add(marks[i], marks[j])) for i in range(len(marks)) for j in range(i + 1, len(marks)))
    h = D0 - md
    # individual primitivity: gcd(z, conj z) (common factor with conjugate <-> divisible by an F_2[t]-prime power pattern)
    gzc = [deg(p_gcd(z, p_conj(z))) for z in marks]
    terms = [sum(1 for c in z if c) for z in marks]
    print(f"{label}: L={L} B={B} raw deg={D0} gcd deg={dg} D={D} max diff(raw)={md} h={h} h/D={h/D:.4f}"
          f" (target {slope}); deg gcd(z,conj z)={gzc}; monomials per mark={terms}")
    return D, h

# her construction: x_j = t^{n_j} (j<=4), x_5 = 1 + sum; H_j spanned by the other four
def her(ns, L):
    xs = [1 << nj for nj in ns]; x5 = 1
    for x in xs: x5 ^= x
    xs.append(x5)
    hyps = []
    for j in range(5):
        H = []
        for c in itertools.product([0, 1], repeat=5):
            if any(c) and c[j] == 0: H.append(c)
        hyps.append(H)
    B = max(ns)
    D, h = analyse(build(xs, L, hyps), f"VII.6 n={ns}", L, B, "8/15=0.5333")
    print(f"      her bound 8L-16B = {8*L-16*B}, referee bound 8L-8B = {8*L-8*B}, D<=15L={15*L}")
for ns, L in (((1, 2, 3, 5), 6), ((1, 2, 3, 5), 20), ((2, 3, 5, 7), 40), ((1, 4, 6, 9), 60)):
    her(ns, L)

# referee's dimension-4 variant: V = span(1, t^a, t^b, t^c); five hyperplanes avoiding 1
def dim4(ns, L, sel):
    xs = [1] + [1 << nj for nj in ns]
    hyps = []
    for lv in sel:
        H = [c for c in itertools.product([0, 1], repeat=4) if any(c) and
             (c[0] ^ (c[1] & lv[0]) ^ (c[2] & lv[1]) ^ (c[3] & lv[2])) == 0]
        hyps.append(H)
    B = max(ns)
    D, h = analyse(build(xs, L, hyps), f"dim4 n={ns}", L, B, "4/7=0.5714")
    print(f"      bound 4L-4B = {4*L-4*B}, D<=7L={7*L}")
sel = [(0,0,0),(1,0,0),(0,1,0),(0,0,1),(1,1,1)]
for ns, L in (((1, 2, 3), 6), ((1, 2, 3), 30), ((2, 5, 7), 80), ((3, 4, 11), 200)):
    dim4(ns, L, sel)
