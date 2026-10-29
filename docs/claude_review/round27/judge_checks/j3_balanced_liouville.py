# Judge check J3 (exact arithmetic): balanced Liouville exponent on the Losev-Manin tower vs dimension.
# k points, 2^{k-1}-1 cuts, each pair separated by 2^{k-2} cuts.  At exact balance n_A = e^w:
#   H_{-K} = R''^2 = e^{(2^{k-1}-1) w};  |t_ab| >= 1  =>  dist(y,e) >> e^{-d_ab/2} = e^{-2^{k-3} w}.
# So the balanced Liouville constant is alpha_bal(k) = (2^{k-1}-1)/2^{k-3}; Vojta predicts >= dim = k-1.
# Cross-check with the record's kappa = (k-1)(1+2 delta)/4 - 1 at the Liouville delta_L(k) = 1/(2^k - 2):
#   A-point exponent gamma = 4/(1+2delta) for every k (H = R''^2, dist <= R''^{-(1/2+delta)}); kappa > 0 iff gamma < dim.
from fractions import Fraction as F
for k in range(3, 10):
    ab = F(2**(k-1)-1, 2**(k-3)); dim = k-1
    dL = F(1, 2**k - 2)
    gamma_L = F(4) / (1 + 2*dL)
    kappa = F(dim)*(1+2*dL)/4 - 1
    print(k, 'alpha_bal =', ab, 'dim =', dim, 'alpha_bal<dim:', ab < dim,
          '| 4/(1+2 delta_L) =', gamma_L, 'equal:', gamma_L == ab, '| kappa(delta_L) =', kappa)
