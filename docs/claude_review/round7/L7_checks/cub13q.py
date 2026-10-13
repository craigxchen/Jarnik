"""EXACT construction of twisted cubics of class cub13 in the Kapranov psi_7 model of Mbar_{0,7}.

Coordinates: P^4 with p_1..p_5 = e_0..e_4, p_6 = (1,1,1,1,1).  Equivalently (integral chart with
the label 7 at infinity) y_6 = 0 and y_{k+1} = x_k: the curve is s -> (y_1(s),...,y_6(s)) and the
Kapranov plane Pi_T (T a 3-subset of [6]) is the triple cluster T^c of the y's.

cub13 = twisted cubic meeting the 13 planes 123,124,125,126,134,135,136,145,246,256,346,356,456,
i.e. the 13 y-triple clusters 456,356,346,345,256,246,236,245,135,134,125,124,123.

Construction (projection from p_1, then from p_2):
  c2 = A (s-1)(s-s45),  c3 = B s(s-1)(s-s35),  c4 = G s(s-s34),   c2(tau)=c3(tau)=c4(tau)=1
  c1 = K (s-s34)(s-s35)(s-s45)  with c1(s36) = c3(s36)
  c0 = cubic interpolating c0 = c2 at s'246, s'256 and c0 = c1 at s'346, s'356
  last condition (plane 456 / cluster 123): exists sigma with c1(sigma)=c2(sigma)!=..., c0(sigma)=c1(sigma).
"""
from fractions import Fraction as Fr
import qpoly as Q


def build(tau, s34, s35, s45):
    tau, s34, s35, s45 = map(Fr, (tau, s34, s35, s45))
    A = 1 / ((tau - 1) * (tau - s45))
    B = 1 / (tau * (tau - 1) * (tau - s35))
    G = 1 / (tau * (tau - s34))
    c2 = Q.scal(A, Q.mul(Q.lin(1), Q.lin(s45)))
    c3 = Q.scal(B, Q.prod([Q.lin(0), Q.lin(1), Q.lin(s35)]))
    c4 = Q.scal(G, Q.mul(Q.lin(0), Q.lin(s34)))
    s36 = (B * s35 + G * s34) / (B * tau)
    K = Q.ev(c3, s36) / ((s36 - s34) * (s36 - s35) * (s36 - s45))
    c1 = Q.scal(K, Q.prod([Q.lin(s34), Q.lin(s35), Q.lin(s45)]))
    s246 = A * s45 / ((A - G) * tau)
    s256 = (A * s45 / B) / tau
    s346 = s35 * s45 / s36
    s356 = K * s34 * s45 / ((K - B) * s36)
    xs = [s246, s256, s346, s356]
    ys = [Q.ev(c2, s246), Q.ev(c2, s256), Q.ev(c1, s346), Q.ev(c1, s356)]
    c0 = Q.interp(xs, ys)
    # quadratic for the 456-plane (cluster {1,2,3} of y): (c1 - c2)/(s - s45)
    R, rem = Q.divmod_(Q.sub(c1, c2), Q.lin(s45))
    assert not rem
    T = Q.sub(c0, c1)
    d = dict(tau=tau, s34=s34, s35=s35, s45=s45, s36=s36, s246=s246, s256=s256, s346=s346,
             s356=s356, A=A, B=B, G=G, K=K, c=[c0, c1, c2, c3, c4], R=R, T=T)
    return d


def E(tau, s34, s35, s45):
    d = build(tau, s34, s35, s45)
    return Q.resultant(d['R'], d['T'])
