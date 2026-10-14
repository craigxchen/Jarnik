"""Index-3 sections P = {y in [-j,j]^b : sum y = 0 mod 3}: alpha = s - rho(Q) versus rho(P).
Exact ranks for (2j+1)^b <= 3200 (exact fmpz), illustrating alpha < rho+1 (contribution far below rho)."""
import itertools, sys
sys.path.insert(0, '.')
from lfi import hilbert
from parity_sections import boxH
for b, j in [(3, 2), (4, 2), (5, 2), (3, 3), (4, 3)]:
    P = [y for y in itertools.product(range(-j, j + 1), repeat=b) if sum(y) % 3 == 0]
    Q = [y for y in itertools.product(range(-j, j + 1), repeat=b) if sum(y) % 3 != 0]
    exact = (2 * j + 1) ** b <= 3200
    HP = hilbert(P, 10 ** 6, cap=2 * j, modp=not exact)
    HQ = hilbert(Q, 10 ** 6, cap=2 * j, modp=not exact)
    HB = boxH(b, j); s = len(HB) - 1
    rP, rQ = len(HP) - 1, len(HQ) - 1
    link = all((HP[e] if e < len(HP) else 0) + (HQ[s - e] if 0 <= s - e < len(HQ) else 0) == HB[e] for e in range(s + 1))
    first = min(e for e in range(s + 1) if (HP[e] if e < len(HP) else 0) < HB[e])
    print(f'mod3 b={b} j={j} |P|={len(P)} s={s} rho(P)={rP} rho(Q)={rQ} alpha=s-rho(Q)={s-rQ} first e with H_P<H_B: {first} linkage={link} {"exact" if exact else "modp"}')
    print('   H_P', HP); sys.stdout.flush()
