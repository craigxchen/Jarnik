# 2-adic classes of odd-norm Gaussian integers; fibre sizes; collision gains.
# Exact integer arithmetic.
from fractions import Fraction
import math

def classes_mod_pow(k):
    # residues of Z[i] modulo (1+i)^k: (1+i)^k Z[i] has index 2^k.
    # represent ideal (1+i)^k: for k even = 2^(k/2) Z[i]; for k odd = 2^((k-1)/2)(1+i) Z[i]
    h = k // 2
    if k % 2 == 0:
        M = 2 ** h
        reps = [(x, y) for x in range(M) for y in range(M)]
        def red(x, y):
            return (x % M, y % M)
    else:
        M = 2 ** h
        # z mod M(1+i): lattice M(1+i)Z[i] = {M(a-b) + i M(a+b)}: x,y both multiples of M with x/M = y/M mod 2
        def red(x, y):
            # reduce y mod 2M? Use canonical: (x mod 2M ... ) simpler: canonical form via lattice basis
            # basis of lattice: (M, M), (M, -M) ; equivalently (2M,0),(M,M)
            # reduce y mod M using (M,M): subtract t*(M,M) with t = y//M
            t = y // M
            x2, y2 = x - t * M, y - t * M
            return (x2 % (2 * M), y2)
        reps = None
    return red

def unit_classes(k):
    red = classes_mod_pow(k)
    S = set()
    B = 2 ** (k // 2 + 2)
    for x in range(-B, B):
        for y in range(-B, B):
            if (x * x + y * y) % 2 == 1:
                S.add(red(x, y))
    return S, red

def fibre_size(k, Nmod):
    # classes mod (1+i)^k of z with N(z) = n for some n = Nmod mod 2^{t}, t large; brute force:
    # count classes c such that some lift z has N(z) == target mod 2^(k+3) for target with target = Nmod mod 8...
    pass

if __name__ == "__main__":
    for k in range(1, 9):
        S, red = unit_classes(k)
        print("k", k, "unit classes mod (1+i)^k:", len(S))
    # fibre: classes realised by points of a fixed circle family; use residues n mod 2^(k+2)
    print("fibre sizes f_k(n): number of unit classes mod (1+i)^k containing z with N(z) = n mod 2^(k+3)")
    for k in range(1, 9):
        red = classes_mod_pow(k)
        T = 2 ** (k + 3)
        B = 2 ** (k + 2)
        fib = {}
        for x in range(B):
            for y in range(B):
                n = x * x + y * y
                if n % 2 == 1:
                    fib.setdefault(n % T, set()).add(red(x, y))
        sizes = sorted(set(len(v) for v in fib.values()))
        print(k, "fibre sizes over n mod 2^(k+3):", sizes)
