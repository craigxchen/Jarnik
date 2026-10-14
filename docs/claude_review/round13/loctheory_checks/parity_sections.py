"""Hilbert functions of lattice sections of the box [-j,j]^b:
  parity circuits P = {y : sum y even}, complement Q = odd points,
and the Hadamard sections Y_t(j) (all 2^t - 1 columns of a t-dim XOR block).
Reports H_P, H_Q, rho(P), rho(Q), alpha = s - rho(Q) (first degree of non-grid generators), and checks the
linkage identity H_P(e) + H_Q(s-e) = H_box(e).
Ranks: exact (fmpz) when EXACT=1 for the instance, else mod p (then rho values are upper bounds, and
alpha = s - rho_p(Q) is a LOWER bound for the true alpha... see write-up §2)."""
import itertools, sys
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from lfi import hilbert

def boxH(b, j):
    h = [1]
    for _ in range(b):
        nh = [0] * (len(h) + 2 * j)
        for a, c in enumerate(h):
            for k in range(2 * j + 1):
                nh[a + k] += c
        h = nh
    return h

def parity(b, j):
    P, Q = [], []
    for y in itertools.product(range(-j, j + 1), repeat=b):
        (P if sum(y) % 2 == 0 else Q).append(y)
    return P, Q

def hadamard_section(t, j):
    """Y_t(j) = {(-k(S_a))_{a != 0} : class sums c in Z^{2^t}, sum c = 0, |k(S_a)| <= j}
    = {y in [-j,j]^{2^t-1} : sum_a y_a chi_a(z) = 0 mod 2^{t-1} for all z}."""
    chars = [a for a in range(1, 2 ** t)]
    zs = list(range(2 ** t))
    par = lambda a, z: bin(a & z).count('1') % 2
    mod = 2 ** (t - 1)
    Y = []
    for y in itertools.product(range(-j, j + 1), repeat=len(chars)):
        ok = True
        for z in zs:
            s = sum(yy * (1 - 2 * par(a, z)) for yy, a in zip(y, chars))
            if s % mod:
                ok = False; break
        if ok:
            Y.append(y)
    return Y

def report(name, P, Q, b, j, modp):
    HB = boxH(b, j)
    s = len(HB) - 1
    HP = hilbert(P, 10 ** 6, cap=2 * j, modp=modp)
    HQ = hilbert(Q, 10 ** 6, cap=2 * j, modp=modp) if Q else [0]
    rP, rQ = len(HP) - 1, len(HQ) - 1
    link = all((HP[e] if e < len(HP) else 0) + (HQ[s - e] if 0 <= s - e < len(HQ) else 0) == HB[e]
               for e in range(s + 1))
    alpha = s - rQ
    print(f'{name}: b={b} j={j} |P|={len(P)} |Q|={len(Q)} s={s} rho(P)={rP} rho(Q)={rQ} alpha=s-rho(Q)={alpha}'
          f' linkage_ok={link} {"modp" if modp else "exact"}')
    print('   H_box', HB)
    print('   H_P  ', HP)
    sys.stdout.flush()

if __name__ == '__main__':
    for b in range(1, 6):
        for j in range(1, 5):
            if (2 * j + 1) ** b > 20000:
                continue
            P, Q = parity(b, j)
            exact = (2 * j + 1) ** b <= 1300
            report('parity', P, Q, b, j, modp=not exact)
    for t, j in [(3, 1), (3, 2)]:
        Y = hadamard_section(t, j)
        b = 2 ** t - 1
        Yset = set(Y)
        Q = [y for y in itertools.product(range(-j, j + 1), repeat=b) if y not in Yset]
        report(f'Y_{t}', Y, Q, b, j, modp=True)
