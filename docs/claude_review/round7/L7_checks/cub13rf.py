"""cub13 construction over Q(x): any of the four parameters may be the variable RF.x()."""
from fractions import Fraction as Fr
from rfun import RF, pscal, pmul, plin, pprod, pev, pdivmod, psub, pinterp, sylvester_res


def build(tau, s34, s35, s45):
    tau, s34, s35, s45 = [a if isinstance(a, RF) else RF(Fr(a)) for a in (tau, s34, s35, s45)]
    A = 1 / ((tau - 1) * (tau - s45))
    B = 1 / (tau * (tau - 1) * (tau - s35))
    G = 1 / (tau * (tau - s34))
    c2 = pscal(A, pmul(plin(RF(1)), plin(s45)))
    c3 = pscal(B, pprod([plin(RF(0)), plin(RF(1)), plin(s35)]))
    c4 = pscal(G, pmul(plin(RF(0)), plin(s34)))
    s36 = (B * s35 + G * s34) / (B * tau)
    K = pev(c3, s36) / ((s36 - s34) * (s36 - s35) * (s36 - s45))
    c1 = pscal(K, pprod([plin(s34), plin(s35), plin(s45)]))
    s246 = A * s45 / ((A - G) * tau)
    s256 = (A * s45 / B) / tau
    s346 = s35 * s45 / s36
    s356 = K * s34 * s45 / ((K - B) * s36)
    xs = [s246, s256, s346, s356]
    ys = [pev(c2, s246), pev(c2, s256), pev(c1, s346), pev(c1, s356)]
    c0 = pinterp(xs, ys)
    R, rem = pdivmod(psub(c1, c2), plin(s45))
    assert not rem
    T = psub(c0, c1)
    return dict(c=[c0, c1, c2, c3, c4], R=R, T=T, s36=s36, K=K, xs=xs)


def E(*args):
    d = build(*args)
    return sylvester_res(d['R'], d['T'])
