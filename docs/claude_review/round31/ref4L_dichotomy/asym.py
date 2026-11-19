# IV.2 family, q fixed, P growing: det E and max delta / (P q^3)
import random, sympy
from math import prod
from fractions import Fraction
from run_tests import iv2_instance
from dich import analyse
x = (-27, -270, 135, 120, 40)
for digits in [50, 100, 200, 400]:
    V, J, P, q = iv2_instance(x, 1, digits)
    o = analyse(V, J, set(), "IV.2", verbose=False)
    print(digits, "detE =", float(o['detE']), " maxdelta/(P q^3) =", float(Fraction(o['maxdelta'], P*q**3)), " signs", o['signs'], "rank", o['rank'])
