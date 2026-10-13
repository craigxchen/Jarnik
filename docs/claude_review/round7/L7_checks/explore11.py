import numpy as np, warnings
warnings.filterwarnings('ignore')
from fractions import Fraction as Fr
import qpoly as Q, param3 as P
from p3curve import genuine
from explore10 import exact, degs, vals
for (name, v), (dn, dd) in degs.items():
    vars_ = ['tau', 'a', 'b', 'c'] if name == 'tabc' else ['mu', 'a', 'b', 'c']
    fn = P.final_tabc if name == 'tabc' else P.final_mabc
    f = lambda x: fn(*[x if w == v else vals[w] for w in vars_])
    num, den = exact(f, dn, dd)
    roots = Pn = np.polynomial.polynomial.polyroots(np.array([float(x) for x in num]))
    gen = []
    for r in roots:
        args = {w: complex(vals[w]) for w in vars_}; args[v] = complex(r)
        if name == 'tabc':
            mu = P.mu_of(args['tau'], args['a'], args['b']); tau = args['tau']
        else:
            mu = args['mu']; tau = P.tau_of(mu, args['a'], args['b'])
        try:
            d = P.data(args['a'], args['b'], args['c'], tau, mu)
        except ZeroDivisionError:
            continue
        if genuine(d): gen.append(r)
    print(name, v, 'genuine', len(gen), np.round(gen, 5), flush=True)
