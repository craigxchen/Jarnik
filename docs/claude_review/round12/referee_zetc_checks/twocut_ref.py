from zonocheck import *
cnt = bad = 0
for d in range(2, 8):
    for a in range(0, d+1):
        for b in range(1, d+1):
            c = d - a - b
            if c < 0: continue
            S1 = tuple(range(0, a+b)); S2 = tuple(range(a, a+b+c))
            for L1 in range(1, a+b):
                for L2 in range(1, b+c):
                    cuts = [(S1, L1), (S2, L2)]
                    A = box(d, cuts)
                    r = rho_p(A, d); z, edges, k0 = zono(d, cuts, A)
                    cnt += 1
                    tag = "OK" if z == r else "GAP"
                    if z != r: bad += 1
                    print(f"a={a} b={b} c={c} L=({L1},{L2}) |A|={len(A)} rho_p={r} zono={z} {tag}", flush=True)
print("blocks", cnt, "gaps", bad)
