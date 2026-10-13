"""degrees (mod p) of the closed-form FINAL in each variable, two parametrisations."""
from modp import M, ratdeg
import param3 as P
vals = dict(tau=M(-2), a=M(3), b=M(5) / 2, c=M(7) / 3, mu=M(11) / 5)
for name, fn, vars_ in [('tabc', P.final_tabc, ['tau', 'a', 'b', 'c']), ('mabc', P.final_mabc, ['mu', 'a', 'b', 'c'])]:
    for v in vars_:
        def f(x):
            args = [x if w == v else vals[w] for w in vars_]
            return fn(*args)
        r = ratdeg(f, 40)
        print(name, v, 'num', r[0], 'den', r[1], flush=True)
