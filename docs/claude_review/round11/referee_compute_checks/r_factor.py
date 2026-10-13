# Factor the basis polynomials P_c(w) = w^{(T..T)} sum_k c_k w^{k'} over Q (flint), classify factors,
# and test whether non-binomial factors are forms in x = w - 1.
import sys, json
import flint
fn = sys.argv[1]
D = json.load(open(fn))
A, others, B = D['A'], D['others'], D['basis']
d = len(others)
ctx = flint.fmpz_mpoly_ctx.get(tuple(f"w{i}" for i in range(d)), 'degrevlex')
xctx = flint.fmpz_mpoly_ctx.get(tuple(f"x{i}" for i in range(d)), 'degrevlex')
gens = ctx.gens()
xg = xctx.gens()
shift = [xg[i] + 1 for i in range(d)]
def poly_of(c):
    terms = {}
    for ck, k in zip(c, A):
        if ck:
            e = tuple(k[x] + 1 for x in others)   # T = 1
            terms[e] = terms.get(e, 0) + ck
    return ctx.from_dict({e: v for e, v in terms.items() if v})
def is_form(p):
    degs = set(sum(e) for e in p.monoms())
    return len(degs) == 1, degs
allinfo = []
for j, c in enumerate(B):
    P = poly_of(c)
    const, facs = P.factor()
    info = {'j': j, 'const': int(const), 'nterms': len(list(P.terms())), 'factors': []}
    order_at_e = 0
    for f, e in facs:
        nt = len(list(f.terms()))
        if nt == 1:
            info['factors'].append(('monomial', str(f), e)); continue
        fx = f.compose(*shift)
        form, degs = is_form(fx)
        mindeg = min(degs)
        order_at_e += mindeg * e
        if nt == 2:
            (m1, c1), (m2, c2) = list(f.terms())
            kind = 'binom w^a - w^b' if c1 == -c2 else f'binom coeffs {c1},{c2}'
            info['factors'].append((kind, str(f), e, 'x-form' if form else 'x-nonform', sorted(degs)))
        else:
            info['factors'].append(('multi', str(fx) if form else str(f), e, 'x-form' if form else 'x-nonform', sorted(degs), nt))
    info['order_at_e'] = order_at_e
    allinfo.append(info)
    print(f"--- basis poly {j}: {info['nterms']} terms, const {info['const']}, order at e = {order_at_e}")
    for fi in info['factors']:
        print("   ", fi)
json.dump(allinfo, open(fn.replace('.json', '_factors.json'), 'w'), default=str)
