from fractions import Fraction as Fr
def F(M): return Fr((M-1)*(M-4),8*(M-2)*(M-3))
def astra_ok(M,D): return (D-1)*F(M) > Fr(3)*D/(2*(M-1))
def rho(M): return Fr(M*M//4, M*(M-1)//2)
def pair_ok(M,D): return D > 4*rho(M)   # pair first moment excludes M-clusters
print("Astra (III.7) exact check:")
for M,D in [(17,4),(18,4),(37,Fr(3,2)),(38,Fr(3,2)),(13,10**9),(14,10**9)]:
    print(M,D,astra_ok(M,D), float((D-1)*F(M)), float(Fr(3)*D/(2*(M-1))))
def first(ok,D):
    for M in range(4,10**5):
        if ok(M,D): return M
    return None
print("D, smallest M excluded: Astra / pair-first-moment")
for D in [Fr(11,10),Fr(5,4),Fr(3,2),Fr(7,4),Fr(2),Fr(21,10),Fr(11,5),Fr(23,10),Fr(5,2),Fr(3),Fr(4),Fr(9,2),Fr(10),Fr(100)]:
    print(float(D), first(astra_ok,D), first(pair_ok,D) if D>2 else None)
# also M=2 for pair: excluded iff D>4
# monotonicity check of astra_ok in M
for D in [Fr(11,10),Fr(3,2),Fr(4)]:
    m0=first(astra_ok,D); assert all(astra_ok(M,D) for M in range(m0,m0+2000)); print("monotone ok",float(D))
