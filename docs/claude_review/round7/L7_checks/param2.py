"""Second parametrisation of the cub13 family (exact). Labels: y_6 = 0, y_2..y_5 cubic forms with
pairwise common roots r_ab ({a,b} in {2,3,4,5}); gauge r45 = 0, r35 = inf, r34 = 1;
b = r24, c = r23, a = r25, tau = cluster 345, mu = cluster 245, lambda_4 = 1.
Stage 2: y_1 - y_2 = kappa (s-p25)(s-p24)(s-sigma), y_1 - y_3 vanishes at p35, p34.
Free parameters here: (b, tau, mu, sp) where sp = sigma' is the root of R = (y2-y3)/(s-c) carrying the
PAIR contact D_23; c is then linear.  Returns the final condition Phi (zero <=> cub13 member)."""
from fractions import Fraction as Fr
import qpoly as Q

def stage1(b, tau, mu, sp):
    a = (tau * mu - b) / (tau + mu - 1 - b)
    l5 = tau + mu - 1 - b
    # c from R(sp) = 0, which is linear in c:
    # (tau - c) mu(mu-1)(sp-a)(sp-b) = tau(tau-b)(sp-1)(mu-a)(mu-c)
    U = mu * (mu - 1) * (sp - a) * (sp - b)
    V = tau * (tau - b) * (sp - 1) * (mu - a)
    c = (tau * U - mu * V) / (U - V)
    l3 = tau * (tau - b) / (tau - c)
    l2 = mu * (mu - 1) / ((mu - a) * (mu - c))
    return dict(a=a, b=b, c=c, tau=tau, mu=mu, sp=sp, l2=l2, l3=l3, l4=Fr(1), l5=l5)

def ys(d):
    a, b, c, l2, l3, l5 = d['a'], d['b'], d['c'], d['l2'], d['l3'], d['l5']
    y5 = Q.scal(l5, Q.prod([Q.lin(0), Q.lin(a)]))
    y4 = Q.prod([Q.lin(0), Q.lin(1), Q.lin(b)])
    y3 = Q.scal(l3, Q.prod([Q.lin(1), Q.lin(c)]))
    y2 = Q.scal(l2, Q.prod([Q.lin(a), Q.lin(b), Q.lin(c)]))
    return y2, y3, y4, y5

def stage2(d):
    a, b, c, tau, mu, sp, l2, l3, l5 = [d[k] for k in ('a', 'b', 'c', 'tau', 'mu', 'sp', 'l2', 'l3', 'l5')]
    p35 = l3 * c / ((l3 - l5) * tau)
    p34 = l3 * c / tau
    p25 = b * c / mu
    p24 = l2 * a * c / ((l2 - 1) * mu)
    # R(s) = l2 s^2 - (l2(a+b) + l3) s + (l2 ab + l3), roots sigma, sp
    sigma = (l2 * a * b + l3) / (l2 * sp)
    d.update(p35=p35, p34=p34, p25=p25, p24=p24, sigma=sigma)
    Phi = (p35 - c) * (p35 - sp) * (p34 - p25) * (p34 - p24) - (p34 - c) * (p34 - sp) * (p35 - p25) * (p35 - p24)
    return Phi

def full(b, tau, mu, sp):
    d = stage1(*map(Fr, (b, tau, mu, sp)))
    Phi = stage2(d)
    return d, Phi

def curve(d):
    """given a member (Phi = 0), return y_1..y_6 as cubic polynomials (forms of degree 3)."""
    y2, y3, y4, y5 = ys(d)
    # y1 - y2 = kappa (s-p25)(s-p24)(s-sigma); kappa from (y1 - y3)(p35) = 0
    base = Q.prod([Q.lin(d['p25']), Q.lin(d['p24']), Q.lin(d['sigma'])])
    kappa = (Q.ev(y3, d['p35']) - Q.ev(y2, d['p35'])) / Q.ev(base, d['p35'])
    y1 = Q.add(y2, Q.scal(kappa, base))
    return [y1, y2, y3, y4, y5, []]
