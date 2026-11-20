# exact-rational recomputation of the union-bound exponent in IX.0 (per unit sJ), worst case J = 40
from fractions import Fraction as F
import math
for name, states in (('integers', 4), ('Gaussian', 7)):
    J = 40
    # s*log(states) <= s*J*log(states)/40 ; coefficient term 9 s log(2P)/(2(b-1)) <= 9 sJ/200 ; saving sJ/5
    coef = math.log(states) / J + 9 / 200 - 1 / 5
    print(name, 'exponent per sJ at J=40:', round(coef, 4), ' margin below -1/20:', round(-1/20 - coef, 4))
# polynomial cutoff: per sJ/A with J >= 40A
print('cutoff', round(math.log(7)/40 + 9/200 - 1/5, 4))
# bounded exponents (twin), per sJ/E, J >= 40 E log(1+3E), b-1 >= 100 E^2 log(2P)/J
print('powers twin', 1/40 + 9/200 - 1/5)
# uniform in C <= exp(sJ/400): extra 9 log C <= 0.0225 sJ
print('uniform-C, integers/Gaussian', round(math.log(4)/40+9/200-1/5+9/400,4), round(math.log(7)/40+9/200-1/5+9/400,4))
