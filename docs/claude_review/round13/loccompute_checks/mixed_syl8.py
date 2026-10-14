from mixed import *
for L, m in ((2, 8), (3, 11)):
    print(f"Syl(8,{L}) m={m}", flush=True)
    print(mixed_b0(*models.sylvester(3, L), m=m, bmax=3))
