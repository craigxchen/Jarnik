"""Equal-weight version of the code test on the double core (one prime of the same degree per column, e=1):
(VI.1) requires minimum Hamming weight of the anchored parity code > (#columns)/4 = q/2.
Exhaustive over all 2^(q-1) codewords (Gray code) for small q. Also first certificate for single Paley q=7 mod 8."""
import sys
sys.path.insert(0,'.')
from profiles import double_core, single_paley, analyse
import io, contextlib
def minwt(S):
    M=len(S); r=len(S[0])
    bits=[[(1+s)//2 for s in row] for row in S]
    rows=[sum((bits[i][j]^bits[0][j])<<j for j in range(r)) for i in range(1,M)]
    best=r+1; w=0; k=len(rows)
    for g in range(1,1<<k):
        i=(g & -g).bit_length()-1
        w^=rows[i]
        if w:
            c=w.bit_count()
            if c<best: best=c
    return best, r
for q in (7,11,19,23):
    b,r=minwt(double_core(q,False))
    print(f"DCo({q}): columns={r}, min Hamming weight={b}, needed > {r/4}: {'PASS' if b>r/4 else 'EXCLUDED (equal weights)'}")
primes=[q for q in range(7, 2000, 8) if all(q % p for p in range(2, int(q ** 0.5) + 1))]
for q in primes:
    buf=io.StringIO()
    with contextlib.redirect_stdout(buf):
        rk,n,cert=analyse(f"P1({q})", single_paley(q))
    if cert is not None:
        print("single Paley q=7 mod 8, first certificate:", buf.getvalue().strip()); break
