# Exact cube tail of one slab statistic: Pr_T(|k(P)| > j), |P| = (q+1)/2, k uniform on {-1,0,1}^q (exact integers).
# Gives mu(T \ A) >= tail (one pattern) and <= 2q * tail (union bound).
import flint, mpmath as mp
mp.mp.dps = 30
out = []
for q, j in [(3251, 257), (5479, 342)]:
    h = (q + 1) // 2
    co = (flint.fmpz_poly([1, 1, 1]) ** h).coeffs()
    tail = sum(int(co[s + h]) for s in range(j + 1, h + 1)) * 2
    lt = mp.log(tail) - h * mp.log(3)
    line = f"q={q} j={j} |P|={h}: ln Pr(|k(P)|>j) = {mp.nstr(lt,10)} ; ln(2q * that) = {mp.nstr(lt+mp.log(2*q),10)}"
    print(line, flush=True); out.append(line)
open('slab_tail_out.txt', 'w').write('\n'.join(out) + '\n')
