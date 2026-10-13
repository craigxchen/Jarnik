# Independent saturation check (referee): V(H_1..H_s) subset union of hyperplanes {x_i=0},{x_i=x_j}?
# Method: Bayer's criterion with python-flint Groebner bases (degrevlex, saturating variable last),
# after a unimodular linear change of coordinates making the linear form the last variable.
import sys, json, time, itertools
import flint
d = 7
X = flint.fmpz_mpoly_ctx.get(tuple(f"x{i}" for i in range(d)), 'degrevlex')
xs = X.gens()
def parse_forms(fn):
    info = json.load(open(fn))
    forms = []   # one entry per usable basis polynomial: list of (form string, exponent)
    for inf in info:
        multis = [f for f in inf['factors'] if f[0] == 'multi']
        if multis and all(f[3] == 'x-form' for f in multis):
            forms.append([(f[1], f[2]) for f in multis])
    return forms
def to_poly(s):
    loc = {f"x{i}": xs[i] for i in range(d)}
    return eval(s.replace('^', '**'), {}, loc)
def change(l):
    """return (phi, phiinv): phi maps x-polys to y-polys (y = same ctx gens) with l(x) -> y_{d-1}"""
    i, j = l
    others = [k for k in range(d) if k != i]
    sigma = {k: idx for idx, k in enumerate(others)}   # others -> 0..d-2
    # phi: x_k -> y_sigma(k) (k != i);  x_i -> y_{d-1} + (y_sigma(j) if j is not None else 0)
    phi = [None] * d
    for k in others:
        phi[k] = xs[sigma[k]]
    phi[i] = xs[d - 1] + (xs[sigma[j]] if j is not None else 0)
    # inverse: y_sigma(k) -> x_k ; y_{d-1} -> x_i - x_j
    inv = [None] * d
    for k in others:
        inv[sigma[k]] = xs[k]
    inv[d - 1] = xs[i] - (xs[j] if j is not None else 0)
    return phi, inv
def strip_last(g):
    # divide by highest power of x_{d-1}
    mons = g.monoms()
    a = min(e[d - 1] for e in mons)
    if a == 0:
        return g
    return X.from_dict({tuple(list(e[:d - 1]) + [e[d - 1] - a]): c for e, c in zip(mons, g.coeffs())})
def saturate(G, l):
    phi, inv = change(l)
    Gy = [g.compose(*phi) for g in G]
    V = flint.fmpz_mpoly_vec(Gy, X)
    GB = V.buchberger_naive()
    assert GB.is_groebner()
    Gs = [strip_last(g) for g in GB]
    Vs = flint.fmpz_mpoly_vec(Gs, X)
    assert Vs.is_groebner()
    Gs = list(Vs.autoreduction())
    return [g.compose(*inv) for g in Gs if not g.is_zero()]
if __name__ == '__main__':
    forms = []
    for lst in parse_forms(sys.argv[1]):
        g = X.from_dict({tuple([0]*d): 1})
        for fs, e in lst:
            g = g * to_poly(fs) ** e
        forms.append(g)   # G_j = product of ALL non-binomial factors of P_j (Lemma 5.1)
    print("forms:", len(forms), [f.total_degree() for f in forms], flush=True)
    G = list(flint.fmpz_mpoly_vec(forms, X).buchberger_naive().autoreduction())
    print("GB of the forms: ", len(G), "elements", flush=True)
    lins = [(i, None) for i in range(d)] + [(i, j) for i, j in itertools.combinations(range(d), 2)]
    t0 = time.time()
    for l in lins:
        G = saturate(G, l)
        unit = any(g.is_constant() and not g.is_zero() for g in G)
        name = f"x{l[0]}" if l[1] is None else f"x{l[0]}-x{l[1]}"
        print(f"after saturating by {name}: {len(G)} GB elements, unit={unit} ({time.time()-t0:.1f}s)", flush=True)
        if unit:
            print("UNIT IDEAL reached: V(forms) is contained in the union of the hyperplanes used so far.")
            break
