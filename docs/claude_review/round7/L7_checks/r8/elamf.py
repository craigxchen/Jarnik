# E_Lambda . f_j for forgetful fibres (label j moves, others fixed), and lambda_A(f_j) = -K.f_j - 2 E.f_j = 2 - 2E.f_j
import random
from fractions import Fraction as Fr
from elamp import degE, md, rnd
from star import LAB
random.seed(11)
for j in range(7):
    out = []
    for rep in range(2):
        mob = []
        for a in range(7):
            if a == j: mob.append([md(rnd()), md(rnd()), 0, 1])
            else: mob.append([0, md(rnd()), 0, 1])
        out.append(degE(mob)[0])
    print('f_' + LAB[j], 'E.f =', out, ' lambda_A(f) =', [2 - 2*d for d in out], flush=True)
