"""ref2_tri: Vieta descent for (A+B-C)^2 = 4AB-4U, positive integers."""
def sol(A,B,C,U): return (A+B-C)**2 == 4*A*B-4*U
def descend(t, U):
    steps = 0
    while True:
        A,B,C = sorted(t)
        assert sol(A,B,C,U)
        if C <= A+B: return (A,B,C), steps
        Cp = 2*(A+B)-C; assert Cp > 0 and Cp*C == (A-B)**2+4*U and A+B+Cp < A+B+C
        t = (A,B,Cp); steps += 1
bad = 0; nseeds = {}
for U in range(1, 201):
    seeds = []
    lim = 4*U  # search box larger than the claimed bound
    for A in range(1, lim+1):
        for B in range(A, lim+1):
            disc = 4*A*B-4*U
            if disc < 0: continue
            r = int(round(disc**0.5))
            if r*r != disc: continue
            for C in {A+B-r, A+B+r}:
                if C >= B and C <= A+B and C > 0:
                    seeds.append((A,B,C))
                    if not (3*A*B <= 4*U and A+B <= 8*U/3): bad += 1
    nseeds[U] = len(set(seeds))
print('seed-bound violations:', bad)
print('number of seeds for U=1..12:', [nseeds[U] for U in range(1,13)])
# ascend from seeds and descend back, e.g. U=4
import itertools
U = 4; cnt = 0
for A in range(1,20):
    for B in range(A,40):
        for C in range(1,200):
            if sol(A,B,C,U):
                s,_ = descend((A,B,C),U); cnt += 1
                A2,B2,C2 = s; assert 3*A2*B2 <= 4*U
print('U=4 solutions with A<20,B<40,C<200 descended:', cnt)
# Astra's U=(uvw)^2 for the computed triangle
print(descend((5668,20381,4553), 4))
