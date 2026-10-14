"""Single full-support congruences P = {y in [-j,j]^b : sum y = 0 mod M}: alpha (first non-grid degree) vs 2jb/M.
alpha computed as first e with H_P(e) < H_box(e); exact ranks when (2j+1)^b <= 3200, else mod p (then the reported
H_P is a lower bound degreewise only cumulatively; first-deviation degree is then a LOWER bound? no: see note)."""
import itertools, sys
sys.path.insert(0, '.')
from lfi import hilbert
from parity_sections import boxH
for M in (4, 8):
    for b, j in [(3, 2), (4, 2), (5, 2), (3, 3), (4, 3), (3, 4), (4, 4), (3, 5), (2, 5), (2, 7)]:
        if (2 * j + 1) ** b > 3200:
            continue
        P = [y for y in itertools.product(range(-j, j + 1), repeat=b) if sum(y) % M == 0]
        HP = hilbert(P, 10 ** 6, cap=2 * j, modp=False)
        HB = boxH(b, j)
        first = min(e for e in range(len(HB)) if (HP[e] if e < len(HP) else 0) < HB[e])
        print(f'M={M} b={b} j={j} |P|={len(P)} rho(P)={len(HP)-1} alpha={first}  2jb/M={2*j*b/M:.2f}  H_P={HP}')
        sys.stdout.flush()
