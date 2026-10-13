"""Judge check J4 (evidence for the tangent-cone form ZE-TC on a T = 1 model): is the real/complex
tangent cone TC(A,m) of compute.md's model M8 (Sylvester H_8, T = 1, L = 2, m = 8; 42 octic forms
Phi_c, c in V(A,8)) contained in a finite union of rational hyperplanes through 0?

Input: compute_checks/had8_T1_L2_m8_forms.json (the author's exact octics; rows over 'mons').
Method: exact Groebner bases over Q (python-flint, degrevlex), Bayer saturation by linear forms
l = sum_i a_i x_i with some a_i = +-1, after the unimodular change of variables making l the last
variable. Linear forms used: x_i (7), x_i - x_j (21), and the tangent forms of the seven
non-collision two-term factors of the author's/referee's M8 bases (u0-u1-u6, u0+u3-u4, u0-u2+u5,
u0-u1+u6, u1+u3-u5, u2+u3-u6, u0-u2-u5 ... see LINS below). If the chain ends in the unit ideal,
V(Phi_c : c) is contained in the union of the hyperplanes used: every point of TC is degenerate.
"""
import json, sys, time, itertools
import flint

HERE = sys.argv[1] if len(sys.argv) > 1 else "../compute_checks/had8_T1_L2_m8_forms.json"
d = 7
X = flint.fmpz_mpoly_ctx.get(tuple(f"x{i}" for i in range(d)), 'degrevlex')
xs = X.gens()


def load():
    info = json.load(open(HERE))
    mons = [tuple(m) for m in info['mons']]
    polys = []
    for row in info['rows']:
        dct = {mons[i]: int(c) for i, c in enumerate(row) if int(c) != 0}
        polys.append(X.from_dict(dct))
    return polys


def change(coeffs):
    """coeffs: list of d ints, some entry +-1. Return (phi, inv): phi substitutes x_i so that
    l(x) becomes y_{d-1}; inv maps back."""
    i = next(k for k in range(d) if abs(coeffs[k]) == 1)
    s = coeffs[i]
    others = [k for k in range(d) if k != i]
    sigma = {k: idx for idx, k in enumerate(others)}
    # l = s x_i + sum_{k != i} a_k x_k  = y_last  =>  x_i = s (y_last - sum a_k y_sigma(k))
    phi = [None] * d
    for k in others:
        phi[k] = xs[sigma[k]]
    rest = X.from_dict({tuple([0] * d): 0})
    for k in others:
        if coeffs[k]:
            rest = rest + coeffs[k] * xs[sigma[k]]
    phi[i] = s * (xs[d - 1] - rest)
    inv = [None] * d
    for k in others:
        inv[sigma[k]] = xs[k]
    lin = X.from_dict({tuple([0] * d): 0})
    for k in range(d):
        if coeffs[k]:
            lin = lin + coeffs[k] * xs[k]
    inv[d - 1] = lin
    return phi, inv


def strip_last(g):
    mons = g.monoms()
    a = min(e[d - 1] for e in mons)
    if a == 0:
        return g
    return X.from_dict({tuple(list(e[:d - 1]) + [e[d - 1] - a]): c for e, c in zip(mons, g.coeffs())})


def saturate(G, coeffs):
    phi, inv = change(coeffs)
    Gy = [g.compose(*phi) for g in G]
    GB = flint.fmpz_mpoly_vec(Gy, X).buchberger_naive()
    Gs = [strip_last(g) for g in GB]
    Gs = list(flint.fmpz_mpoly_vec(Gs, X).buchberger_naive().autoreduction())
    return [g.compose(*inv) for g in Gs if not g.is_zero()]


def vec(*pairs):
    v = [0] * d
    for k, c in pairs:
        v[k] += c
    return v


LINS = [("x%d" % i, vec((i, 1))) for i in range(d)]
LINS += [("x%d-x%d" % (i, j), vec((i, 1), (j, -1))) for i, j in itertools.combinations(range(d), 2)]
LINS += [
    ("x0-x1-x6", vec((0, 1), (1, -1), (6, -1))),   # w0 - w1 w6
    ("x0+x6-x1", vec((0, 1), (6, 1), (1, -1))),    # w0 w6 - w1
    ("x0-x2-x5", vec((0, 1), (2, -1), (5, -1))),   # w0 - w2 w5
    ("x0+x5-x2", vec((0, 1), (5, 1), (2, -1))),    # w0 w5 - w2
    ("x0+x3-x4", vec((0, 1), (3, 1), (4, -1))),    # w0 w3 - w4
    ("x1+x3-x5", vec((1, 1), (3, 1), (5, -1))),    # w1 w3 - w5
    ("x2+x3-x6", vec((2, 1), (3, 1), (6, -1))),    # w2 w3 - w6
    ("x1+x6-x0", vec((1, 1), (6, 1), (0, -1))),    # w1 w6 - w0 (referee basis)
    ("x2+x5-x0", vec((2, 1), (5, 1), (0, -1))),    # w2 w5 - w0 (referee basis)
]

if __name__ == "__main__":
    t0 = time.time()
    F = load()
    print("forms:", len(F), "degrees", sorted(set(f.total_degree() for f in F)), flush=True)
    G = list(flint.fmpz_mpoly_vec(F, X).buchberger_naive().autoreduction())
    print(f"GB of the 42 octics: {len(G)} elements ({time.time()-t0:.1f}s)", flush=True)
    for name, cf in LINS:
        G = saturate(G, cf)
        unit = any(g.is_constant() and not g.is_zero() for g in G)
        print(f"after saturating by {name}: {len(G)} GB elements, unit={unit} ({time.time()-t0:.1f}s)", flush=True)
        if unit:
            print("UNIT IDEAL reached: TC(M8) is contained in the union of the hyperplanes used so far.")
            break
    else:
        print("chain ended WITHOUT the unit ideal; remaining GB degrees:", sorted(g.total_degree() for g in G)[:20])
