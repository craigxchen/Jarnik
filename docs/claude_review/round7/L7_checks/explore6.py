from modp import M, ratdeg, roots_in_set
import param2 as P2
base = dict(b=M(3), tau=M(-2), mu=M(5) / 2, sp=M(7) / 3)
def full(b, tau, mu, sp):
    d = P2.stage1(b, tau, mu, sp); return P2.stage2(d)
for var in ['sp', 'b', 'tau', 'mu']:
    def f(x):
        a = dict(base); a[var] = x
        return full(a['b'], a['tau'], a['mu'], a['sp'])
    r = ratdeg(f, 40)
    print(var, 'deg num', r[0], 'deg den', r[1], flush=True)
