from rlib import *
import itertools
# slices Lambda_j: rho <= s-|j| ?
for s in range(1, 8):
    out = []
    for j in range(0, s+1):
        Lam = [k for k in cube_points(s) if sum(k) == j]
        E = Exact(Lam)
        out.append((j, len(Lam), E.rho, s-j))
        assert E.rho <= s - j
    print("s", s, [(j, n, r, b) for j, n, r, b in out], flush=True)
# author's indicator degree: prod_Z (1-k^2) prod_N (k^2-k)/2 has degree 2|Z|+2|N| (> s-j when Z nonempty)
# corrected indicator prod_Z (1-k) prod_N (k^2-k)/2: check it is the indicator of kappa on Lambda_j
for s in range(1, 7):
    for j in range(0, s+1):
        Lam = [k for k in cube_points(s) if sum(k) == j]
        for kap in Lam:
            Z = [x for x in range(s) if kap[x] == 0]; N = [x for x in range(s) if kap[x] == -1]
            for k in Lam:
                v = 1
                for x in Z: v *= (1 - k[x])
                for x in N: v *= (k[x]**2 - k[x]) // 2
                assert v == (1 if k == kap else 0)
            assert len(Z) + 2*len(N) == s - j
print("corrected indicator verified s<=6")
# Cor 2.2(3): S={0,1}, L=1, d=6
A = slabbox(6, [((0,1), 1)])
E = Exact(A)
print("Cor2.2: |A|", len(A), "rho", E.rho, "H", E.H)
for u in [(1,2,3,5,0,0), (2,7,-3,11,0,0), (1,2,1,1,0,0), (1,-1,3,5,0,0), (1,2,3,5,7,0)]:
    print(u, "nu", E.nu(u))
# width of k(S) and widths of coordinate directions
print("width k(S):", max(k[0]+k[1] for k in A) - min(k[0]+k[1] for k in A))
