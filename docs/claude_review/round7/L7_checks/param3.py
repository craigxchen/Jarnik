"""cub13 family, closed form (exact; works over any field: Fraction, modp.M, RF, complex).

Gauge (y-chart, y_6 = 0, label 7 at infinity): clusters 456 at s=0, 356 at s=inf, 346 at s=1.
Parameters: a = r_256, b = r_246, c = r_236 (common roots of y_2 with y_5, y_4, y_3),
tau = r_345, mu = r_245, related by  a (tau + mu - 1 - b) = tau mu - b.
   y5 = l5 s(s-a),  y4 = s(s-1)(s-b),  y3 = l3 (s-1)(s-c),  y2 = l2 (s-a)(s-b)(s-c),  y6 = 0,
   l5 = tau+mu-1-b,  l3 = tau(tau-b)/(tau-c),  l2 = mu(mu-1)/((mu-a)(mu-c)).
Label-1 clusters: 135 at p35, 134 at p34, 125 at p25, 124 at p24, 123 at sigma.
The pair contact of 2 and 3 is at sp = iota(c), iota = the involution of P^1 exchanging p35<->p34 and
p25<->p24.  Membership condition:  R(sp) = 0,  R = (y2-y3)/(s-c) = l2(s-a)(s-b) - l3(s-1).
FINAL(...) below is (numerator-free) l2 (sp-a)(sp-b) - l3 (sp-1)."""


def mu_of(tau, a, b):
    return (b + a * (tau - 1 - b)) / (tau - a)


def tau_of(mu, a, b):
    return (b + a * (mu - 1 - b)) / (mu - a)


def a_of(tau, mu, b):
    return (tau * mu - b) / (tau + mu - 1 - b)


def data(a, b, c, tau, mu):
    l5 = tau + mu - 1 - b
    l3 = tau * (tau - b) / (tau - c)
    l2 = mu * (mu - 1) / ((mu - a) * (mu - c))
    p34 = c * (tau - b) / (tau - c)
    p35 = (tau - b) * c / (tau * (tau - b) - l5 * (tau - c))
    p25 = b * c / mu
    p24 = (mu - 1) * a * c / (mu * (mu - 1) - (mu - a) * (mu - c))
    S1, P1 = p35 + p34, p35 * p34
    S2, P2 = p25 + p24, p25 * p24
    sp = (S1 * P2 - S2 * P1 - c * (P2 - P1)) / ((P2 - P1) + c * (S1 - S2))
    # sigma = other root of R(s) = l2 s^2 - (l2(a+b) + l3) s + (l2 a b + l3)
    sigma = (l2 * a * b + l3) / (l2 * sp)
    final = l2 * (sp - a) * (sp - b) - l3 * (sp - 1)
    return dict(a=a, b=b, c=c, tau=tau, mu=mu, l2=l2, l3=l3, l5=l5, p34=p34, p35=p35, p25=p25,
                p24=p24, sp=sp, sigma=sigma, final=final)


def final_tabc(tau, a, b, c):
    return data(a, b, c, tau, mu_of(tau, a, b))['final']


def final_mabc(mu, a, b, c):
    return data(a, b, c, tau_of(mu, a, b), mu)['final']
