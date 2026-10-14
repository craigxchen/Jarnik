from zonocheck import box, rho_p, zono
from tc_enum import blocks
cnt = gaps = 0
for dd, cols, Ls in blocks(6):
    cuts = list(zip(cols, Ls)); A = box(dd, cuts)
    r = rho_p(A, dd); z, edges, k0 = zono(dd, cuts, A); cnt += 1
    if z != r: gaps += 1; print("GAP", dd, cols, Ls, "rho_p", r, "zono", z, flush=True)
    if cnt % 50 == 0: print("done", cnt, "gaps", gaps, flush=True)
print("three-cut blocks", cnt, "gaps", gaps, flush=True)
# remaining two-cut type (a,b,c) = (3,1,3) at d = 7
for L1 in (1, 2, 3):
    for L2 in (1, 2, 3):
        cuts = [(tuple(range(0, 4)), L1), (tuple(range(3, 7)), L2)]
        A = box(7, cuts); r = rho_p(A, 7); z, e, k0 = zono(7, cuts, A)
        print(f"twocut a=3 b=1 c=3 L=({L1},{L2}) rho_p={r} zono={z} {'OK' if r == z else 'GAP'}", flush=True)
